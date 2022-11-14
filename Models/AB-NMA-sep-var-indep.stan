// AB-NMA model with covariance selection
// horseshoe prior on treatment random effects standard deviations

functions{
  matrix make_Gamma(vector gamma, int ntrt){
    matrix[ntrt, ntrt] Gamma;
    
  for(i in 2:ntrt){
    Gamma[i, i] = 1;
    for(j in 1:(i - 1)){
      Gamma[i, j] = gamma[(i - 1) * (i - 2) / 2 + j];
      Gamma[j, i] = 0;
    }
  }
  Gamma[1, 1] = 1;
    
    return Gamma;
  }
}

data {
  int<lower = 0> len; // total number of data size
  int<lower = 0> ntrt ; // total number of treatments
  int<lower = 0> nstudy ; // total number of studies
  int<lower = 0> t[len]; // indices of treatment
  int<lower = 0> s[len]; // indices of study
  int<lower = 0> r[len]; // number of events
  int<lower = 0> totaln [len]; // number of participants
  int<lower = 0, upper = 1> higher_better; // to indicate whether higher absolute risk is better or not
  vector<lower = 0>[ntrt] zeros; // mean effects in MVN
  
  real<lower = 0> scale_global_lambda;  // scale for half-t on tau for lambdas
  
  real<lower = 1> nu_global_lambda; // df for half-t on tau for lambda
  real<lower = 1> nu_local_lambda; // df for half-t on omega for lambdas
  
  real<lower = 0> slab_scale_lambda; // slab scale for regularized horseshoe
  real<lower = 0> slab_df_lambda; // slab df for regularized horseshoe
}

transformed data{
  int n_gamma = (ntrt * (ntrt - 1)) / 2; // number of free elements in Gamma
}
parameters {
  vector[nstudy] beta; // random effects for study
  matrix [nstudy, ntrt] vi; // random study-treatment effects
  vector [ntrt] mu; // fixed effect for treatment
  real<lower = 0> sigma_beta; // variance of study-level random effect
  
  real<lower = 0> tau_lambda;  // global shrinkage parameter for lambda params
  vector<lower = 0>[ntrt] omega_lambda;  // local shrinkage parameter for lambda params
  real<lower = 0> c2_lambda; 
  
  vector<lower = 0>[ntrt] z_lambda;
  
}
transformed parameters {
  vector<lower = 0, upper = 1>[len] p; // Study_specified Absolute Risk

  real<lower = 0> c_lambda;
  
  vector<lower = 0>[ntrt] omega_lambda_tilde; // regularized 
  // vector<lower = 0>[ntrt] scale_gamma; // adjust scale for gammas based on lambda
  
  vector<lower = 0>[ntrt] lambda;  // elements of Lambda matrix in cholesky decomp
  // vector[ntrt] gamma;  // elements of gamma matrix in cholesky decomp


  for(i in 1:len){
    p[i] = inv_logit(mu[t[i]] + beta[s[i]] + vi[s[i], t[i]]);
  }
  
  c_lambda = slab_scale_lambda * sqrt(c2_lambda);
  
  omega_lambda_tilde = sqrt(c_lambda ^ 2 * square(omega_lambda) ./ (c_lambda ^ 2 + tau_lambda ^ 2 * square(omega_lambda)));
  
  // lambda contains random effects standard deviations
  lambda = z_lambda .* omega_lambda_tilde * tau_lambda;
  
}
model {
  // Regularized horseshoe parameterization for elements of Lambda, Gamma
  z_lambda ~ normal(0, 1);
  
  tau_lambda ~ student_t(nu_global_lambda, 0, scale_global_lambda);
  omega_lambda ~ student_t(nu_local_lambda, 0, 1);
  
  c2_lambda ~ inv_gamma(0.5 * slab_df_lambda, 0.5 * slab_df_lambda);
  
  sigma_beta ~ cauchy(0, 0.5);

  for(j in 1:ntrt){
    mu[j] ~ normal(0, 10);
  }
  for(j in 1:nstudy){
    beta[j] ~ normal(0, sigma_beta);
    for(i in 1:ntrt){
      vi[j, i] ~ normal(0, lambda[i]);
    }
  }
    for(i in 1:len){
    r[i] ~ binomial (totaln[i], p[i]);
  }
}
generated quantities {
  real AR[ntrt]; // Populaiton averaged absolute risk
  real MU[ntrt];
  real SD[ntrt]; // standard deviation of covariance matrix
  real OR[ntrt, ntrt]; // Odds Ratio
  real LOR[ntrt , ntrt]; // Log odds Ratio
  real cLOR [ntrt, ntrt]; // cLog odds Ratio
  // real CORR [ntrt, ntrt]; // correlation matrix
  // cov_matrix[ntrt] Sigma = multiply_lower_tri_self_transpose(Sigma_chol);
  int rk[ntrt]; // rank of AR
  int rank_prob [ntrt, ntrt]; // rank probability

  real sumll = 0; // sum of log likelihood
  real log_likelihood [len];
  real likelihood [len];
  for(j in 1:ntrt){
    AR[j] = inv_logit(mu[j]/ sqrt(1 + (sigma_beta ^ 2 + lambda[j] ^ 2) * 256 / 75 / pi() / pi()));
    SD[j] = lambda[j];
    MU[j] = mu[j];
  }
  for(j in 1:ntrt){
    for(k in 1:ntrt){
      OR[j,k] = AR[j] / (1 - AR[j]) / AR[k] * (1 - AR[k]);
      LOR[j,k] = log(OR[j,k]);
      cLOR [j,k] = MU[j] - MU[k];
      // CORR [j,k] = Sigma[j,k] / SD[j] / SD[k];
    }
  }
  for(j in 1:ntrt){
    rk[j]= higher_better?ntrt-rank(AR,j):rank(AR,j)+1;
  }
  for(j in 1:ntrt){
    for(k in 1:ntrt){
      rank_prob[k,j]=(rk[k]==j);
    }
  }
  for(i in 1:len){
    log_likelihood[i]=r[i] * log(p[i]) + (totaln[i] - r[i]) * log(1 - p[i]);
    likelihood[i] = exp(log_likelihood[i]);
    sumll = sumll + log_likelihood[i];
  }
}

