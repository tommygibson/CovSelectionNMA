// sample from prior cov matrix with regularized horseshoe in NMA

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
  int ntrt;
  
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
  int n_gamma = (ntrt * (ntrt - 1)) / 2;
}

parameters {
  real<lower = 0> tau_lambda;  // global shrinkage parameter for lambda params
  real<lower = 0> tau_gamma;  // global shrinkage parameter for gamma params
  vector<lower = 0>[ntrt] omega_lambda;  // local shrinkage parameter for lambda params
  vector<lower = 0>[n_gamma] omega_gamma;   // local shrinkage parameter for gamma params
  real<lower = 0> c2_lambda; 
  real<lower = 0> c2_gamma;
  
  vector<lower = 0>[ntrt] z_lambda;
  vector[n_gamma] z_gamma;
}

transformed parameters {
  real<lower = 0> c_lambda;
  real<lower = 0> c_gamma;
  
  vector<lower = 0>[ntrt] omega_lambda_tilde; // regularized 
  vector<lower = 0>[n_gamma] omega_gamma_tilde; // regularized
  
  vector<lower = 0>[ntrt] lambda;  // elements of Lambda matrix in cholesky decomp
  vector[n_gamma] gamma;  // elements of gamma matrix in cholesky decomp
  
  matrix[ntrt, ntrt] Gamma;
  cov_matrix[ntrt] Sigma;
  
  c_lambda = slab_scale_lambda * sqrt(c2_lambda);
  c_gamma = slab_scale_gamma * sqrt(c2_gamma);
  
  omega_lambda_tilde = sqrt(c_lambda ^ 2 * square(omega_lambda) ./ (c_lambda ^ 2 + tau_lambda ^ 2 * square(omega_lambda)));
  omega_gamma_tilde = sqrt(c_gamma ^ 2 * square(omega_gamma) ./ (c_gamma ^ 2 + tau_gamma ^ 2 * square(omega_gamma)));
  
  lambda = z_lambda .* omega_lambda_tilde * tau_lambda;
  gamma = z_gamma .* omega_gamma_tilde * tau_gamma;
  Gamma = make_Gamma(gamma, ntrt);
  Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply(lambda, Gamma));
  
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

}

generated quantities {
  vector[ntrt] SD;
  matrix[ntrt, ntrt] CORR;
  
  for(i in 1:ntrt){
    SD[i] = sqrt(Sigma[i, i]);
  }
  for(i in 1:ntrt){
    for(j in 1:ntrt){
      CORR[i, j] = Sigma[i, j] / SD[i] / SD[j];
    }
  }
}

