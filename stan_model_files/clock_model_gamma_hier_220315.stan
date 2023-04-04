data
{
  int<lower=1> N; // number of observations
  vector<lower=0>[N] y;    // genetic distances for pairwise sequences of pair ij
  vector<lower=0>[N] x;  // time elapsed for pair ij
  int<lower=1> K; // number of pairs
  array[N] int<lower=1,upper=K> pair_idx; // pair id for each observation
  int<lower=0> N_pr; // number of distances to predict
  vector<lower=0>[N_pr] x_pr; // time elapsed for distances to predict
}
transformed data
{
  real log_evolutionary_rate_pol_within_host_avg = log(10^(-2.5));
  real log_evolutionary_rate_pol_within_host_sd = 0.2;
}
parameters
{
  real log_alpha1;
  vector[K] log_alpha1_pair;
  real<lower=0> log_alpha1_pair_sd;
  real log_phi;
  vector[K] log_phi_pair;
  real<lower=0> log_phi_pair_sd;
}
transformed parameters
{
  vector[N] mu;
  vector[N] alpha;
  vector[N] beta;

  beta = exp(-(log_phi + log_phi_pair[pair_idx]));
  mu = exp(log_alpha1 + log_alpha1_pair[pair_idx] ) .* x;
  alpha = mu .* beta;
}
model
{
  // priors
  target += normal_lpdf(log_alpha1 | log_evolutionary_rate_pol_within_host_avg, log_evolutionary_rate_pol_within_host_sd);
  target += normal_lpdf(log_alpha1_pair | 0, log_alpha1_pair_sd);
  target += exponential_lpdf(log_alpha1_pair_sd | 10);
  target += normal_lpdf(log_phi | 0, 5);
  target += normal_lpdf(log_phi_pair | 0, log_phi_pair_sd);
  target += exponential_lpdf(log_phi_pair_sd | 10);

  // likelihood
  target += gamma_lpdf(y | alpha , beta);
}
generated quantities
{
  array[N] real log_lik;
  vector[N_pr] mu_pr;
  real<lower=0> beta_pr;
  array[N_pr] real y_pr;

  // pointwise log lik
  for(i in 1:N)
  {
    log_lik[i] = gamma_lpdf(y[i] | alpha[i] , beta);
  }

  // distances predicted for particular x_pr
  mu_pr = exp(log_alpha1) * x_pr;
  beta_pr = exp(-log_phi);
  for(i in 1:N_pr)
  {
    y_pr[i] = gamma_rng( mu_pr[i]*beta_pr, beta_pr );
  }
}
