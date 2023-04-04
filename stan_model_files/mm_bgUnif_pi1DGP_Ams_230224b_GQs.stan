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
  int<lower=0> Nk_max; //Max no .observ's per age
  int<lower=0> Nj_max; //Max no .observ's per recip 1yr age
  int<lower=0> Nj_max_coarse; //Max no .observ's per recip age band
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

  // data summaries
  real uxsd_src = sd(ux_src);
  real uxsd_rec = sd(ux_rec);

  // HSGP basis functions for f at unique inputs
  matrix[Nux_src, M[1]] hsgp_eigenvalues_src = PHI(N, M[1], L[1], ux_src);
  matrix[Nux_rec, M[1]] hsgp_eigenvalues_rec = PHI(N, M[1], L[1], ux_rec);
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
  vector[M[1]] hspg_spectral_dens_rec;
  vector[N] hsgp_f_src;
  vector[N] hsgp_f_rec;

  // spectral densities for f
  hspg_spectral_dens_src = diagSPD_EQ(gpscale, lscale, L[1], M[1]);
  hspg_spectral_dens_rec = diagSPD_EQ(gpscale, lscale, L[2], M[2]);

  // HSGP approximation to f at inputs
  hsgp_f_src = (hsgp_eigenvalues_src * (hspg_spectral_dens_src .* hsgp_beta))[x_2_ux_src];
  hsgp_f_rec = (hsgp_eigenvalues_rec * (hspg_spectral_dens_rec .* hsgp_beta))[x_2_ux_rec];

  // mixing weights
  y_mix = inv_logit(logit_y_mix_0 + hsgp_f_src + hsgp_f_rec);

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
