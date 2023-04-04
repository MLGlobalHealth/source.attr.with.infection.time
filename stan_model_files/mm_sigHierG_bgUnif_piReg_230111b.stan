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

model
{
  target += normal_lpdf(log_alpha1_pair | 0, log_alpha1_pair_sd);
  target += normal_lpdf(log_phi_pair | 0, log_phi_pair_sd);
  target += normal_lpdf(logit_y_mix_0 | 0, 2); //-5.293305, 5
  target += normal_lpdf(beta_age_src | 0, 1);
  target += normal_lpdf(beta_age_rcp | 0, 1);


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
