// 3RE model WITHOUT covariance selection
// Horseshoe prior on elements of special cholesky decomposition

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
  vector<lower = 0>[3] lambda;  // elements of Lambda matrix in cholesky decomp
  vector[3] gamma;  // elements of gamma matrix in cholesky decomp
  real<lower = 0> tau_lambda;  // global shrinkage parameter for lambda params
  real<lower = 0> tau_gamma;  // global shrinkage parameter for gamma params
  vector<lower = 0>[3] omega_lambda;  // local shrinkage parameter for lambda params
  vector<lower = 0>[3] omega_gamma;   // local shrinkage parameter for gamma params

}

transformed parameters {
  vector<lower = 0, upper = 1>[S] pi1;
  vector<lower = 0, upper = 1>[S] pi0;
  vector<lower = 0, upper = 1>[S] psi;
  matrix[3, 3] Gamma;
  cholesky_factor_cov[3] Sigma_chol;
  
  for(i in 1:S){
    pi1[i] = inv_logit(theta[i, 1] + theta[i, 2] / 2);
    pi0[i] = inv_logit(theta[i, 1] - theta[i, 2] / 2);
    psi[i] = inv_logit(theta[i, 3]);
  }
  Gamma = make_Gamma(gamma);
  Sigma_chol = diag_pre_multiply(lambda, Gamma);
}

model {
  
  // hyperpriors
  theta0[1] ~ normal(a, b);
  theta0[2] ~ normal(c, d);
  theta0[3] ~ normal(f, g);
  
  tau_lambda ~ cauchy(0, 1);
  tau_gamma ~ cauchy(0, 1);
  omega_lambda ~ cauchy(0, 1);
  omega_gamma ~ cauchy(0, 1);
  
  lambda ~ normal(0, tau_lambda * omega_lambda);
  gamma ~ normal(0, tau_gamma * omega_gamma);

  for(i in 1:S){
    theta[i,] ~ multi_normal_cholesky(theta0, Sigma_chol);
  }
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

