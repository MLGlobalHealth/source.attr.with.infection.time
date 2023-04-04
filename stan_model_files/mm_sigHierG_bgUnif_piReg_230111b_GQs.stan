data
{
  int<lower=1> N; //number of observations
  int<lower=1> P; //number of unique recipients
  int<lower=1> S; //number of source categories
  vector<lower=0>[N] y; // genetic distances for pairwise sequences of pair ij
  vector<lower=0>[N] x; // time elapsed for pair ij
  int<lower=0> A; // number of binned age groups
  int<lower=1,upper=A> obs_to_age_src_idx[N];
  int<lower=1,upper=A> obs_to_age_rcp_idx[N];
  int<lower=1> pt_idx[N]; // unique recipient ID for pair ij
  int pt_map[P,N]; // idx of ID for pair ij
  matrix[S,N] idx_src; // idx of source categories
  matrix[S,N] idx_rec; // idx of recipient categories
  real log_alpha1;
  real log_alpha1_pair_sd;
  real log_phi;
  real log_phi_pair_sd;
}

transformed data
{
  real<lower=0> unif_pdf = 1/max(y);
  real unif_lpdf = log(unif_pdf);
}

parameters
{
  vector[N] log_alpha1_pair;
  vector[N] log_phi_pair;
  real logit_y_mix_0;
  vector[A] beta_age_src; // age coefficient of mixture probability
  vector[A] beta_age_rcp; // age coefficient of mixture probability
}

transformed parameters
{
  vector<lower=0, upper=1>[N] y_mix; // probability that observation is a transmission pair
  vector[N] tpair_mean;
  vector[N] tpair_alpha;
  vector[N] tpair_beta;

  //y_mix = inv_logit(logit_y_mix_0 + beta_age[obs_to_age_idx]);
  y_mix = inv_logit(logit_y_mix_0 + beta_age_src[obs_to_age_src_idx] + beta_age_rcp[obs_to_age_rcp_idx]);
  tpair_beta = exp(-(log_phi + log_phi_pair));
  tpair_mean = exp(log_alpha1 + log_alpha1_pair ) .* x;
  tpair_alpha = tpair_mean .* tpair_beta;
}

generated quantities
{
  vector[N] log_lik;

  for (i in 1:N){
      log_lik[i] = log_mix(y_mix[i], gamma_lpdf(y[i] | tpair_alpha[i], tpair_beta[i]), unif_lpdf);
  }

}
