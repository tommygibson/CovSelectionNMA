// 3RE model with covariance selection
// Regularized Horseshoe prior on elements of special cholesky decomposition
// Covariance matrix for random effects is decomposed as
// Sigma = Lambda * Gamma * t(Gamma) * Lambda
// Lambda = diag(lambda[1], lambda[2], lambda[3])
// Gamma = lower triangular, 1 on diagonal, gamma[1,2], gamma[1,3], gamma[2,3]
// so 3 lambda parameters and 3 gamma parameters, each with separate 
// regularized horseshoe priors

functions{
  matrix make_Gamma(vector gamma){
    int K = num_elements(gamma);
    matrix[3, 3] Gamma;
    
    for(i in 1:3){
      Gamma[i, i] = 1;
    }
    Gamma[2, 1] = gamma[1];
    Gamma[3, 1] = gamma[2];
    Gamma[3, 2] = gamma[3];
    Gamma[1, 2] = 0;
    Gamma[1, 3] = 0;
    Gamma[2, 3] = 0;
    
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
  real<lower = 0> scale_global_gamma; // scale for half-t on tau for gammas
  
  real<lower = 1> nu_global_lambda; // df for half-t on tau for lambda
  real<lower = 1> nu_global_gamma; // df for half-t on tau for gammas
  real<lower = 1> nu_local_lambda; // df for half-t on omega for lambdas
  real<lower = 1> nu_local_gamma; // df for half-t on omega for gammas
  
  real<lower = 0> slab_scale_lambda; // slab scale for regularized horseshoe
  real<lower = 0> slab_df_lambda; // slab df for regularized horseshoe
  real<lower = 0> slab_scale_gamma; // slab scale for regularized horseshoe
  real<lower = 0> slab_df_gamma; // slab df for regularized horseshoe
  
}

transformed data {
  int N[S];
  for(i in 1:S){
    N[i] = n1[i] + n0[i];
  }
}

parameters {
  matrix[S, 3] theta; // random effects for log odds event, log(OR), log odds RF
  vector[3] theta0; // mean hyperparameter vector
  real<lower = 0> tau_lambda;  // global shrinkage parameter for lambda params
  real<lower = 0> tau_gamma;  // global shrinkage parameter for gamma params
  vector<lower = 0>[3] omega_lambda;  // local shrinkage parameter for lambda params
  vector<lower = 0>[3] omega_gamma;   // local shrinkage parameter for gamma params
  real<lower = 0> c2_lambda; 
  real<lower = 0> c2_gamma;
  
  vector<lower = 0>[3] z_lambda;
  vector[3] z_gamma;

}

transformed parameters {
  vector<lower = 0, upper = 1>[S] pi1;
  vector<lower = 0, upper = 1>[S] pi0;
  vector<lower = 0, upper = 1>[S] psi;
  
  real<lower = 0> c_lambda;
  real<lower = 0> c_gamma;
  
  vector<lower = 0>[3] omega_lambda_tilde; // regularized 
  vector<lower = 0>[3] omega_gamma_tilde; // regularized
  
  vector<lower = 0>[3] lambda;  // elements of Lambda matrix in cholesky decomp
  vector[3] gamma;  // elements of gamma matrix in cholesky decomp
  
  matrix[3, 3] Gamma;
  cholesky_factor_cov[3] Sigma_chol;
  
  for(i in 1:S){
    pi1[i] = inv_logit(theta[i, 1] + theta[i, 2] / 2);
    pi0[i] = inv_logit(theta[i, 1] - theta[i, 2] / 2);
    psi[i] = inv_logit(theta[i, 3]);
  }
  c_lambda = slab_scale_lambda * sqrt(c2_lambda);
  c_gamma = slab_scale_gamma * sqrt(c2_gamma);
  
  omega_lambda_tilde = sqrt(c_lambda ^ 2 * square(omega_lambda) ./ (c_lambda ^ 2 + tau_lambda ^ 2 * square(omega_lambda)));
  omega_gamma_tilde = sqrt(c_gamma ^ 2 * square(omega_gamma) ./ (c_gamma ^ 2 + tau_gamma ^ 2 * square(omega_gamma)));
  
  lambda = z_lambda .* omega_lambda_tilde * tau_lambda;
  gamma = z_gamma .* omega_gamma_tilde * tau_gamma;
  Gamma = make_Gamma(gamma);
  Sigma_chol = diag_pre_multiply(lambda, Gamma);
}

model {
  
  // Regularized horseshoe parameterization for elements of Lambda, Gamma
  z_lambda ~ normal(0, 1);
  z_gamma ~ normal(0, 1);
  
  tau_lambda ~ student_t(nu_global_lambda, 0, scale_global_lambda);
  tau_gamma ~ student_t(nu_global_gamma, 0, scale_global_gamma);
  omega_lambda ~ student_t(nu_local_lambda, 0, 1);
  omega_gamma ~ student_t(nu_local_gamma, 0, 1);
  
  c2_lambda ~ inv_gamma(0.5 * slab_df_lambda, 0.5 * slab_df_lambda);
  c2_gamma ~ inv_gamma(0.5 * slab_df_gamma, 0.5 * slab_df_gamma);

  for(i in 1:S){
    theta[i,] ~ multi_normal_cholesky(theta0, Sigma_chol);
  }
  
  // mean hyperpriors
  theta0[1] ~ normal(a, b);
  theta0[2] ~ normal(c, d);
  theta0[3] ~ normal(f, g);
  
  // binomial likelihood
  y1 ~ binomial(n1, pi1);
  y0 ~ binomial(n0, pi0);
  n1 ~ binomial(N, psi);
  
}

generated quantities {
  vector[3] SD;
  matrix[3, 3] CORR;
  matrix[3, 3] Sigma = multiply_lower_tri_self_transpose(Sigma_chol);
  
  for(i in 1:3){
    SD[i] = sqrt(Sigma[i, i]);
  }
  for(i in 1:3){
    for(j in 1:3){
      CORR[i, j] = Sigma[i, j] / SD[i] / SD[j];
    }
  }
}

