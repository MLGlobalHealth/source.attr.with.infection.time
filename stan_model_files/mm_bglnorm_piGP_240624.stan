functions{
  vector diagSPD_EQ(real alpha, real rho, real L, int M)
    {
      return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
    }

  // Eigenfunctions
  matrix PHI(int M, real L, vector x)
  {
    return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
  }

  matrix kron_mvprod(matrix A, matrix B, matrix V)
  {
    return (B*V) * transpose(A);
  }

matrix hsgp(int A, real alpha, real rho1, real rho2, real L1, real L2, int M1, int M2,
            matrix PHI1, matrix PHI2, matrix z)
{
  vector[M1] sqrt_spd_1 = diagSPD_EQ(alpha, rho1, L1, M1);
  vector[M2] sqrt_spd_2 = diagSPD_EQ(alpha, rho2, L2, M2);

  matrix[A,A] f = kron_mvprod(
    diag_post_multiply( PHI1, sqrt_spd_1 ),
    diag_post_multiply( PHI2, sqrt_spd_2 ),
    z
  );

  return(f);
}

}

data{
  int<lower=1> N; //No. observations or no.pairs
  int<lower=1> P; //number of unique recipients
  int<lower=1> S; //number of source categories
  vector<lower=0>[N] y; //Genetic Distances
  vector<lower=0>[N] x; //Time elapsed

  int<lower=1> D; //Number of Dimensions
  vector[D] L; //Boundaries for each dimension
  array[D] int M;  //Number of basis functions in each dimension
  int<lower=1> M_nD; //M1*M2

  matrix[M_nD, D] indices; //Matrix of D-tuples
  int<lower=0> A; // number of age groups
  int<lower=0> n; // number of 1yr ages
  matrix[n, D] ages; //Age matrix

  // parameters from clock model
  real log_alpha1;
  real log_alpha1_pair_sd;
  real log_phi;
  real log_phi_pair_sd;

  //Indexing
  array[N, D] int<lower=0,upper=N> coordinates; // indexing for ages of sources and recipients
  array[N] int<lower=1> pt_idx; // unique recipient ID for pair ij
  array[P,N] int pt_map; // idx of ID for pair ij

}

transformed data{
  matrix[n, M[1]] PHI1;
  matrix[n, M[2]] PHI2;
  PHI1 = PHI(M[1], L[1], ages[,1]);
  PHI2 = PHI(M[2], L[2], ages[,2]);
}

parameters{
  vector[N] log_alpha1_pair;
  vector[N] log_phi_pair;
  real logit_y_mix_0; //Mixing probability parameter (on logit scale)

  vector<lower=0>[D] lscale; //Lengthscale
  real<lower=0> gpscale; //Amplitude
  matrix[M[1],M[2]] z1;
  
  real epsilon_base; // baseline coefficient of mean of background mixture component
  real epsilon; // slope coefficient of mean of background mixture component
  real<lower=0> sigma; // standard deviation of background mixture component
}

transformed parameters{
  matrix[n,n] f; //Latent function
  {
        f = hsgp(n, gpscale, lscale[1], lscale[2],
              L[1], L[2], M[1], M[2], PHI1, PHI2, z1);
  }
  vector<lower=0, upper=1>[N] y_mix; //Mixing probability parameter (on natural scale)
  vector[N] tpair_mean; //Mean of transmission pair
  vector[N] tpair_alpha; //Shape param of signal component
  vector[N] tpair_beta; //Rate param of signal component
  vector[N] mu_background; // mean of the background

  for(i in 1:N){
      y_mix[i] = inv_logit(logit_y_mix_0 + f[coordinates[i,1],coordinates[i,2]]);
  }
  tpair_beta = exp(-(log_phi + log_phi_pair)); //Beta = e^(-phi0_log + phipair_log) where phipair is a vector i.e log(beta) = -log(phi0) + log(phipair)
  tpair_mean = exp(log_alpha1 + log_alpha1_pair) .* x; //mu = e^(log_alpha1 + log_alpha_pair) so log(mu) = log_alpha1 + log_alpha_pair
  tpair_alpha = tpair_mean .* tpair_beta; //For gamma dist, mu = alpha / beta -> alpha = mu .* beta

  mu_background = epsilon_base + epsilon*x;

}

model{
  target += normal_lpdf(log_alpha1_pair | 0, log_alpha1_pair_sd);
  target += normal_lpdf(log_phi_pair | 0, log_phi_pair_sd);
  target += normal_lpdf(logit_y_mix_0 | 0, 2); //mean = 0.005 = 200/200^2, (was 0 for 0.5) -5.293305, 5

  target += normal_lpdf(epsilon_base | -2, 0.5);
  target += normal_lpdf(epsilon | 0, 0.2);
  target += exponential_lpdf(sigma | 2);

  // GP priors
  lscale[1] ~ inv_gamma(5, 5);
  lscale[2] ~ inv_gamma(5, 5);
  gpscale ~ normal(0,0.15);

 for(i in 1:M[1]){
    for(j in 1:M[2]){
      z1[i,j] ~ normal(0,1);
    }
  }

  for (i in 1:N){
    target += log_mix(y_mix[i], gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]), lognormal_lpdf(y[i] | mu_background[i],sigma)); //Background pairs given uniform
  }
}

generated quantities{
  int npairs;
  int k;
  real num;     // create vector of appropriate size for all npairs
  vector[2] den_one;
  vector<lower=0,upper=1>[N] tpair_prob;
  vector<lower=0,upper=1>[N] tpair_prob_w;

  for (i in 1:N)
  {

    npairs = sum(pt_map[pt_idx[i],]); //count number of pairs for recipient in pair i

    vector[npairs+1] den;     // create vector of appropriate size for all npairs

    // numerator of pr(transm pair)
    num = gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]) +
    log(y_mix[i]); // first add the llk of pair i being a true pair

    k = 1;
    if(npairs>1){
      for(j in 1:N)
      { // add llk of other pairs for same recipient being non-pairs
          if(i!=j && pt_map[pt_idx[i],j]==1){
            num += lognormal_lpdf(y[i] | mu_background[i],sigma) + log1m(y_mix[i]);
            k += 1;
          }
       }
    }

    // denominator of pr(transm pair)
    k = 1;
    if(npairs>1){
      for(j in 1:N)
      { // all combinations of pairs with other non-pairs
          if(pt_map[pt_idx[i],j]==1){
              den[k] = gamma_lpdf(y[j] | tpair_alpha[j], tpair_beta[j]) + log(y_mix[i]);
              den[npairs+1] = lognormal_lpdf(y[i] | mu_background[i],sigma) + log1m(y_mix[i]);
            for(l in 1:N)
            {
              if(j!=l && pt_map[pt_idx[i],l]==1){
                den[k] += lognormal_lpdf(y[i] | mu_background[i],sigma) + log1m(y_mix[i]);
                den[npairs+1] += lognormal_lpdf(y[i] | mu_background[i],sigma) + log1m(y_mix[i]);
              }
           }
            k += 1;
          }
       }
    }else{ // if only one possible pair for the recipient, only add pr(pair | non-pair)
      den_one[1] = gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]) + log(y_mix[i]);
      den_one[2] = lognormal_lpdf(y[i] | mu_background[i],sigma) + log1m(y_mix[i]);
    }

    if(npairs==1){
        tpair_prob_w[i] = exp(num-log_sum_exp(den_one));
    }else{
        tpair_prob_w[i] = exp(num-log_sum_exp(den));
    }

    tpair_prob[i] = exp(
          gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]) +
          log(y_mix[i]) -
          (
            log_mix( y_mix[i],
              gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]),
              lognormal_lpdf(y[i] | mu_background[i],sigma)
              )
          )
        );

  }
}
