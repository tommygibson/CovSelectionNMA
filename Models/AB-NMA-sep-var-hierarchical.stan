// AB-NMA model with covariance selection
// horseshoe prior on elements of special cholesky decomposition

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
  real<lower = 0> scale_gamma;
  
  real<lower = 1> nu_global_lambda; // df for half-t on tau for lambda
  real<lower = 1> nu_local_lambda; // df for half-t on omega for lambdas
  
  real<lower = 0> slab_scale_lambda; // slab scale for regularized horseshoe
  real<lower = 0> slab_df_lambda; // slab df for regularized horseshoe
}

transformed data{
  int n_gamma = (ntrt * (ntrt - 1)) / 2; // number of free elements in Gamma
}
parameters {
  vector[nstudy] beta_aux; // random effects for study
  vector [ntrt] mu_aux; // fixed effect for treatment
  matrix [nstudy, ntrt] vi; // random study-treatment effects
  real mu0; // mean of treatments
  real<lower = 0> sigma_mu; // sd of treatments
  real<lower = 0> sigma_beta; // variance of study-level random effect
  
  real<lower = 0> tau_lambda;  // global shrinkage parameter for lambda params
  vector<lower = 0>[ntrt] omega_lambda;  // local shrinkage parameter for lambda params
  real<lower = 0> c2_lambda; 
  
  vector<lower = 0>[ntrt] z_lambda;
  vector[n_gamma] z_gamma;
  
}
transformed parameters {
  vector<lower = 0, upper = 1>[len] p; // Study_specified Absolute Risk
  matrix[nstudy, ntrt] vi0; // mean of study-treatment random effects
  vector[ntrt] mu;
  vector[nstudy] beta;

  real<lower = 0> c_lambda;
  
  vector<lower = 0>[ntrt] omega_lambda_tilde; // regularized 
  // vector<lower = 0>[ntrt] scale_gamma; // adjust scale for gammas based on lambda
  
  vector<lower = 0>[ntrt] lambda;  // elements of Lambda matrix in cholesky decomp
  // vector[ntrt] gamma;  // elements of gamma matrix in cholesky decomp
  
  matrix[ntrt, ntrt] Gamma_unscaled; // Unscaled gamma matrix
  matrix[ntrt, ntrt] Gamma_scales;  // scaling factors for gamma matrix conditional on lambdas
  matrix[ntrt, ntrt] Gamma;         // scaled gamma matrix
  cholesky_factor_cov[ntrt] Sigma_chol;
  
  for(i in 1:nstudy){
    beta[i] = sigma_beta * beta_aux[i]; // non-centered parameterization for beta random effects
    for(j in 1:ntrt){
      mu[j] = mu0 + sigma_mu * mu_aux[j]; // non-centered parameterization for mu random effects
      
      vi0[i, j] = mu[j] + beta[i];
    }
  }
  for(i in 1:len){
    p[i] = inv_logit(vi[s[i], t[i]]);
  }
  
  Gamma_scales[ntrt, ntrt] = 1;
  
  c_lambda = slab_scale_lambda * sqrt(c2_lambda);
  
  omega_lambda_tilde = sqrt(c_lambda ^ 2 * square(omega_lambda) ./ (c_lambda ^ 2 + tau_lambda ^ 2 * square(omega_lambda)));
  
  lambda = z_lambda .* omega_lambda_tilde * tau_lambda;
  
  for(i in 1:(ntrt - 1)){
    Gamma_scales[i, i] = 1;
    for(j in (i + 1):ntrt){
      Gamma_scales[j, i] = scale_gamma * sqrt((1 / (1 / lambda[i]^2 + 1 / lambda[j]^2)));
      Gamma_scales[i, j] = 0;
    }
  }
  // gamma = z_gamma .* scale_gamma;
  Gamma_unscaled = make_Gamma(z_gamma, ntrt);
  Gamma = Gamma_unscaled .* Gamma_scales;
  Sigma_chol = diag_pre_multiply(lambda, Gamma);
  
}
model {
  // Regularized horseshoe parameterization for elements of Lambda, Gamma
  z_lambda ~ normal(0, 1);
  z_gamma ~ normal(0, 1);
  
  tau_lambda ~ student_t(nu_global_lambda, 0, scale_global_lambda);
  omega_lambda ~ student_t(nu_local_lambda, 0, 1);
  
  c2_lambda ~ inv_gamma(0.5 * slab_df_lambda, 0.5 * slab_df_lambda);
  
  mu0 ~ normal(0, 10);
  sigma_mu ~ cauchy(0, 1);
  
  sigma_beta ~ cauchy(0, 1);

  mu_aux ~ normal(0, 1);
  beta_aux ~ normal(0, 1);
  
  for(j in 1:nstudy){
    vi[j,]~ multi_normal_cholesky(vi0[j,], Sigma_chol);
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
  real CORR [ntrt, ntrt]; // correlation matrix
  cov_matrix[ntrt] Sigma = multiply_lower_tri_self_transpose(Sigma_chol);
  int rk[ntrt]; // rank of AR
  int rank_prob [ntrt, ntrt]; // rank probability

  real sumll = 0; // sum of log likelihood
  real log_likelihood [len];
  real likelihood [len];
  for(j in 1:ntrt){
    AR[j] = inv_logit(mu[j]/ sqrt(1 + (sigma_beta ^ 2 + Sigma [j,j]) * 256 / 75 / pi() / pi()));
    SD[j] = sqrt(Sigma[j,j]);
    MU[j] = mu[j];
  }
  for(j in 1:ntrt){
    for(k in 1:ntrt){
      OR[j,k] = AR[j] / (1 - AR[j]) / AR[k] * (1 - AR[k]);
      LOR[j,k] = log(OR[j,k]);
      cLOR [j,k] = MU[j] - MU[k];
      CORR [j,k] = Sigma[j,k] / SD[j] / SD[k];
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

