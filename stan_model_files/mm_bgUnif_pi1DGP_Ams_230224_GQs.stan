functions {
  vector diagSPD_EQ(real alpha, real rho, real L, int M) {
    return alpha * sqrt(sqrt(2*pi()) * rho) * exp(-0.25*(rho*pi()/2/L)^2 * linspaced_vector(M, 1, M)^2);
  }
  matrix PHI(int M, real L, vector x) {
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

  int<lower=1> Nux; //number of unique observed ages
  vector[Nux] ux; // unique ages
  array[N] int<lower=1, upper=Nux> x_2_ux; // index of unique values

  int<lower=1> D; //Number of Dimensions
  //vector[D] L; //Boundaries for each dimension
  int<lower=1> M;    // number of basis functions
  real<lower=1> c;

  int<lower=0> A; // number of age groups

  // parameters from clock model
  real log_alpha1;
  real log_alpha1_pair_sd;
  real log_phi;
  real log_phi_pair_sd;

    //Indexing
  int<lower=0> Nk_max; //Max obs per age
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

  // define boundary value
	real L = c*max(ux);
  
  // HSGP basis functions for f at unique inputs
  matrix[Nux, M] hsgp_eigenvalues = PHI(M, L, ux);
}

parameters
{
  // parameters of meta-population clock model
  vector[N] log_alpha1_pair;
  vector[N] log_phi_pair;
  real logit_y_mix_0; // baseline (logit scale)

  // mixing parameters
  vector[M] hsgp_beta; // the HSGP basis functions coefficients
  real<lower=0> lscale; // lengthscale of the squared exp kernel
  real<lower=0> gpscale; // sqrt variance scale of the squared exp kernel
}

transformed parameters
{
  vector<lower=0, upper=1>[N] y_mix; // mixing prob parameter depending on (natural scale)
  vector[N] tpair_mean;
  vector[N] tpair_alpha;
  vector[N] tpair_beta;
  vector[M] hspg_spectral_dens;
  vector[Nux] f;

  // spectral densities for f
  hspg_spectral_dens = diagSPD_EQ(gpscale, lscale, L, M);

  // HSGP approximation to f at inputs
  f = (hsgp_eigenvalues * (hspg_spectral_dens .* hsgp_beta));

  // mixing weights
  y_mix = inv_logit(logit_y_mix_0 + f[x_2_ux]);

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
  vector[N] log_lik;

  for (i in 1:N){
      log_lik[i] = log_mix(y_mix[i], gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]), unif_lpdf);
  }

}
