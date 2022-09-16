// 3RE model with covariance selection
// Regularized Horseshoe prior on elements of special cholesky decomposition
// Covariance matrix for random effects is decomposed as
// Sigma = Lambda * Gamma * t(Gamma) * Lambda
// Lambda = diag(lambda[1], lambda[2], lambda[3])
// Gamma = lower triangular, 1 on diagonal, gamma[1,2], gamma[1,3], gamma[2,3]
// so 3 lambda parameters and 3 gamma parameters, each with separate 
// regularized horseshoe priors

functions{
  matrix make_Gamma(vector gamma, int n_re){
    matrix[n_re, n_re] Gamma;
    
  for(i in 2:n_re){
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
  int<lower=0> S;  // number of studies
  int y1[S];  // number of events in risk factor group
  int y0[S];  // number of events in no risk factor group
  int n1[S];  // number of subjects in risk factor group
  int n0[S];  // number of subjects in non risk factor group
  
  real a;  // prior mean of beta0
  real b; // prior sd of beta0
  real c; // prior mean of delta0
  real d; // prior sd of delta0
  real f; // prior mean of nu0
  real g; // prior sd of nu0
  real<lower = 0> scale_global_lambda;  // scale for half-t on tau for lambdas
  real<lower = 0> scale_gamma;
  
  real<lower = 1> nu_global_lambda; // df for half-t on tau for lambda
  real<lower = 1> nu_local_lambda; // df for half-t on omega for lambdas
  
  real<lower = 0> slab_scale_lambda; // slab scale for regularized horseshoe
  real<lower = 0> slab_df_lambda; // slab df for regularized horseshoe
  
}

transformed data {
  int N[S];
  int n_gamma;
  int n_re = 3;
  n_gamma = (n_re * (n_re - 1)) / 2;
  for(i in 1:S){
    N[i] = n1[i] + n0[i];
  }
  
}

parameters {
  matrix[S, n_re] theta; // random effects for log odds event, log(OR), log odds RF
  vector[n_re] theta0; // mean hyperparameter vector
  real<lower = 0> tau_lambda;  // global shrinkage parameter for lambda params
  vector<lower = 0>[n_re] omega_lambda;  // local shrinkage parameter for lambda params
  real<lower = 0> c2_lambda; 
  
  vector<lower = 0>[n_re] z_lambda;
  vector[n_gamma] z_gamma;

}

transformed parameters {
  vector<lower = 0, upper = 1>[S] pi1;
  vector<lower = 0, upper = 1>[S] pi0;
  vector<lower = 0, upper = 1>[S] psi;
  
  real<lower = 0> c_lambda;
  
  vector<lower = 0>[n_re] omega_lambda_tilde; // regularized 
  // vector<lower = 0>[n_re] scale_gamma; // adjust scale for gammas based on lambda
  
  vector<lower = 0>[n_re] lambda;  // elements of Lambda matrix in cholesky decomp
  // vector[n_re] gamma;  // elements of gamma matrix in cholesky decomp
  vector[n_gamma] gamma; // scaled gamma
  matrix[n_re, n_re] Gamma;         // scaled gamma matrix
  cholesky_factor_cov[n_re] Sigma_chol;
  
  for(i in 1:S){
    pi1[i] = inv_logit(theta[i, 1] + theta[i, 2] / 2);
    pi0[i] = inv_logit(theta[i, 1] - theta[i, 2] / 2);
    psi[i] = inv_logit(theta[i, 3]);
  }
  
  c_lambda = slab_scale_lambda * sqrt(c2_lambda);
  
  omega_lambda_tilde = sqrt(c_lambda ^ 2 * square(omega_lambda) ./ (c_lambda ^ 2 + tau_lambda ^ 2 * square(omega_lambda)));
  
  lambda = z_lambda .* omega_lambda_tilde * tau_lambda;
  
  gamma = scale_gamma * z_gamma; // scaling z_gamma
  Gamma = make_Gamma(gamma, n_re);
  Sigma_chol = diag_pre_multiply(lambda, Gamma);
}

model {
  
  // Regularized horseshoe parameterization for elements of Lambda
  z_lambda ~ normal(0, 1);
  // elements of Gamma
  z_gamma ~ normal(0, 1);
  
  tau_lambda ~ student_t(nu_global_lambda, 0, scale_global_lambda);
  omega_lambda ~ student_t(nu_local_lambda, 0, 1);
  
  c2_lambda ~ inv_gamma(0.5 * slab_df_lambda, 0.5 * slab_df_lambda);

  // mean hyperpriors
  theta0[1] ~ normal(a, b);
  theta0[2] ~ normal(c, d);
  theta0[3] ~ normal(f, g);
  
  for(i in 1:S){
    theta[i,] ~ multi_normal_cholesky(theta0, Sigma_chol);
  }
  
  // binomial likelihood
  y1 ~ binomial(n1, pi1);
  y0 ~ binomial(n0, pi0);
  n1 ~ binomial(N, psi);
  
}

generated quantities {
  vector[n_re] SD;
  matrix[n_re, n_re] CORR;
  matrix[n_re, n_re] Sigma = multiply_lower_tri_self_transpose(Sigma_chol);
  
  for(i in 1:n_re){
    SD[i] = sqrt(Sigma[i, i]);
  }
  for(i in 1:n_re){
    for(j in 1:n_re){
      CORR[i, j] = Sigma[i, j] / SD[i] / SD[j];
    }
  }
}

