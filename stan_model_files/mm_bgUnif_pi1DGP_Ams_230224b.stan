functions {
  vector diagSPD_EQ(real alpha, real rho, real L, int M) {
    return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
  }
  matrix PHI(int N, int M, real L, vector x) {
    return sin(diag_post_multiply(rep_matrix(pi()/(2*L) * (x+L), M), linspaced_vector(M, 1, M)))/sqrt(L);
  }
}

data
{
  int<lower=1> N; //number of observations
  vector<lower=0>[N] y; // genetic distances for pairwise sequences of pair ij
  vector<lower=0>[N] x; // time elapsed for pair ij
  int<lower=1> P; //number of unique recipients
  int<lower=1> S; //number of source categories

  int<lower=1> Nux_src; //number of unique observed ages
  int<lower=1> Nux_rec; //number of unique observed ages
  vector[Nux_src] ux_src; // unique ages
  vector[Nux_rec] ux_rec; // unique ages
  array[N] int<lower=1, upper=Nux_src> x_2_ux_src; // index of unique values
  array[N] int<lower=1, upper=Nux_rec> x_2_ux_rec; // index of unique values

  int<lower=1> D; //Number of Dimensions
  vector[D] L; //Boundaries for each dimension
  int<lower=1> M[D];    // number of basis functions

  int<lower=0> A; // number of age groups

  // parameters from clock model
  real log_alpha1;
  real log_alpha1_pair_sd;
  real log_phi;
  real log_phi_pair_sd;

    //Indexing
  int<lower=0> Nk_max; //Max obs per age
  int<lower=0,upper=N> coordinates[N, D]; // indexing for ages of sources and recipients
  int<lower=0> Nj_max; //Max obs per recip 1yr age
  int<lower=0> Nj_max_coarse; //Max obs per recip age band
  int<lower=0> age_recip_coarse[Nj_max_coarse, A, A];
  int<lower=1> pt_idx[N]; // unique recipient ID for pair ij
  int pt_map[P,N]; // idx of ID for pair ij
  matrix[S,N] idx_src; // idx of source categories
  matrix[S,N] idx_rec; // idx of recipient categories

}

transformed data
{
  real<lower=0> unif_pdf = 1/max(y);
  real unif_lpdf = log(unif_pdf);

  // HSGP basis functions for f at unique inputs
  matrix[Nux_src, M[1]] hsgp_eigenvalues_src = PHI(Nux_src, M[1], L[1], ux_src);
  matrix[Nux_rec, M[2]] hsgp_eigenvalues_rec = PHI(Nux_rec, M[2], L[2], ux_rec);
}

parameters
{
  // parameters of meta-population clock model
  vector[N] log_alpha1_pair;
  vector[N] log_phi_pair;
  real logit_y_mix_0; // baseline (logit scale)

  // mixing parameters
  vector[M[1]] hsgp_beta; // the HSGP basis functions coefficients
  real<lower=0> lscale; // lengthscale of the squared exp kernel
  real<lower=0> gpscale; // sqrt variance scale of the squared exp kernel
}

transformed parameters
{
  vector<lower=0, upper=1>[N] y_mix; // mixing prob parameter depending on (natural scale)
  vector[N] tpair_mean;
  vector[N] tpair_alpha;
  vector[N] tpair_beta;
  vector[M[1]] hspg_spectral_dens_src;
  vector[M[2]] hspg_spectral_dens_rec;
  vector[Nux_src] hsgp_f_src;
  vector[Nux_rec] hsgp_f_rec;

  // spectral densities for f
  hspg_spectral_dens_src = diagSPD_EQ(gpscale, lscale, L[1], M[1]);
  hspg_spectral_dens_rec = diagSPD_EQ(gpscale, lscale, L[2], M[2]);

  // HSGP approximation to f at inputs
  hsgp_f_src = (hsgp_eigenvalues_src * (hspg_spectral_dens_src .* hsgp_beta));
  hsgp_f_rec = (hsgp_eigenvalues_rec * (hspg_spectral_dens_rec .* hsgp_beta));

  // mixing weights
  y_mix = inv_logit(logit_y_mix_0 + hsgp_f_src[x_2_ux_src] + hsgp_f_rec[x_2_ux_rec]);

  // transformed parameters of meta-population clock model
  tpair_beta = exp(-(log_phi + log_phi_pair));
  tpair_mean = exp(log_alpha1 + log_alpha1_pair ) .* x;
  tpair_alpha = tpair_mean .* tpair_beta;
}

model
{
  // priors of meta-population clock model
  target += normal_lpdf(log_alpha1_pair | 0, log_alpha1_pair_sd);
  target += normal_lpdf(log_phi_pair | 0, log_phi_pair_sd);

  // priors of mixing weights
  target += normal_lpdf(logit_y_mix_0 | 0, 2);
  target += normal_lpdf(hsgp_beta | 0, 1);
  target += inv_gamma_lpdf(lscale | 5, 5);
  target += normal_lpdf(gpscale | 0, .5);

  for (i in 1:N)
  {
    target += log_mix( y_mix[i],
                       gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]),
                       unif_lpdf
                       );
  }
}

generated quantities
{
  vector[N] d_pred;
  int npairs;
  int k;
  real num;     // create vector of appropriate size for all npairs
  vector[2] den_one;
  vector<lower=0,upper=1>[N] tpair_prob;
  vector<lower=0,upper=1>[N] tpair_prob_w;
  row_vector<lower=0,upper=1>[N] idx;
  matrix<lower=0>[S,S] flows; //rows=sources
  matrix<lower=0>[S,S] pflows; // prop. total flows between groups
  matrix<lower=0>[S,S] pflows_to; // prop. flows within recipients B from each source
  vector<lower=0>[S] pflows_from; // prop. flows from sources A


  for (i in 1:N)
  {

    // generate distances using parameters for posterior predictive check
    d_pred[i] = log_sum_exp(log(y_mix[i]) + gamma_rng(tpair_alpha[i], tpair_beta[i]), log1m(y_mix[i]) + uniform_rng(0,max(y)));

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
            num += unif_lpdf + log1m(y_mix[i]);
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
              den[npairs+1] = unif_lpdf + log1m(y_mix[i]);
            for(l in 1:N)
            {
              if(j!=l && pt_map[pt_idx[i],l]==1){
                den[k] += unif_lpdf + log1m(y_mix[i]);
                den[npairs+1] += unif_lpdf + log1m(y_mix[i]);
              }
           }
            k += 1;
          }
       }
    }else{ // if only one possible pair for the recipient, only add pr(pair | non-pair)
      den_one[1] = gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]) + log(y_mix[i]);
      den_one[2] = unif_lpdf + log1m(y_mix[i]);
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
              unif_lpdf
              )
          )
        );

  }
  for(a in 1:S){
    for(b in 1:S){
        idx = (idx_src[a,] .* idx_rec[b,]);
        flows[a,b] = idx * tpair_prob_w;
      }
  }
    pflows = flows/sum(flows[:,:]);
  for(a in 1:S){
    for(b in 1:S){
        pflows_to[a,b] = flows[a,b]/sum(flows[:,b]);
    }
  }
  for(a in 1:S){
    pflows_from[a] = sum(flows[a,:])/sum(flows);
  }
}
