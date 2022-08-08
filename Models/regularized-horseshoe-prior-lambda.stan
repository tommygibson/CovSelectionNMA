// sample from prior cov matrix with regularized horseshoe in NMA

functions{
  matrix make_Gamma(vector gamma, vector lambda, int ntrt){
    matrix[ntrt, ntrt] Gamma;
    
  for(i in 2:ntrt){
    Gamma[i, i] = 1;
    for(j in 1:(i - 1)){
      if(lambda[i] < 0.025 || lambda[j] < 0.025){
        Gamma[i, j] = 0.1 * gamma[(i - 1) * (i - 2) / 2 + j];
        Gamma[j, i] = 0;
      } else{
      Gamma[i, j] = gamma[(i - 1) * (i - 2) / 2 + j];
      Gamma[j, i] = 0;
      }
    }
  }
  Gamma[1, 1] = 1;
    
    return Gamma;
  }
}

data {
  int ntrt;
  
  real<lower = 0> scale_global_lambda;  // scale for half-t on tau for lambdas
  real<lower = 0> scale_gamma; // scale for half-t on tau for gammas
  
  real<lower = 1> nu_global_lambda; // df for half-t on tau for lambda
  real<lower = 1> nu_local_lambda; // df for half-t on omega for lambdas
  
  real<lower = 0> slab_scale_lambda; // slab scale for regularized horseshoe
  real<lower = 0> slab_df_lambda; // slab df for regularized horseshoe
}

transformed data {
  int n_gamma = (ntrt * (ntrt - 1)) / 2;
}

parameters {
  real<lower = 0> tau_lambda;  // global shrinkage parameter for lambda params
  vector<lower = 0>[ntrt] omega_lambda;  // local shrinkage parameter for lambda params
  real<lower = 0> c2_lambda; 
  
  vector<lower = 0>[ntrt] z_lambda;
  vector[n_gamma] z_gamma;
}

transformed parameters {
  real<lower = 0> c_lambda;
  
  vector<lower = 0>[ntrt] omega_lambda_tilde; // regularized 
  
  vector<lower = 0>[ntrt] lambda;  // elements of Lambda matrix in cholesky decomp
  vector[n_gamma] gamma;  // elements of gamma matrix in cholesky decomp
  
  matrix[ntrt, ntrt] Gamma;
  cov_matrix[ntrt] Sigma;
  
  c_lambda = slab_scale_lambda * sqrt(c2_lambda);
  
  omega_lambda_tilde = sqrt(c_lambda ^ 2 * square(omega_lambda) ./ (c_lambda ^ 2 + tau_lambda ^ 2 * square(omega_lambda)));
  
  lambda = z_lambda .* omega_lambda_tilde * tau_lambda;
  gamma = z_gamma * scale_gamma;
  Gamma = make_Gamma(gamma, lambda, ntrt);
  Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply(lambda, Gamma));
  
}

model {
  // Regularized horseshoe parameterization for elements of Lambda, Gamma
  z_lambda ~ normal(0, 1);
  z_gamma ~ normal(0, 1);
  
  tau_lambda ~ student_t(nu_global_lambda, 0, scale_global_lambda);
  omega_lambda ~ student_t(nu_local_lambda, 0, 1);
  
  c2_lambda ~ inv_gamma(0.5 * slab_df_lambda, 0.5 * slab_df_lambda);

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

