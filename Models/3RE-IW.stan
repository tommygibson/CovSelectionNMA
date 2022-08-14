// 3RE model WITHOUT covariance selection
// inverse-Wishart prior on covariance matrix



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
  matrix[3, 3] D = diag_matrix(rep_vector(1.0, 3));
  for(i in 1:S){
    N[i] = n1[i] + n0[i];
  }
  
}

parameters {
  matrix[S, 3] theta; // random effects for log odds event, log(OR), log odds RF
  vector[3] theta0; // mean hyperparameter vector
  cov_matrix[3] Sigma;

}

transformed parameters {
  vector<lower = 0, upper = 1>[S] pi1;
  vector<lower = 0, upper = 1>[S] pi0;
  vector<lower = 0, upper = 1>[S] psi;
  
  for(i in 1:S){
    pi1[i] = inv_logit(theta[i, 1] + theta[i, 2] / 2);
    pi0[i] = inv_logit(theta[i, 1] - theta[i, 2] / 2);
    psi[i] = inv_logit(theta[i, 3]);
  }
}

model {
  // binomial likelihood
  y1 ~ binomial(n1, pi1);
  y0 ~ binomial(n0, pi0);
  n1 ~ binomial(N, psi);
  
  // priors

  
  // hyperpriors
  theta0[1] ~ normal(a, b);
  theta0[2] ~ normal(c, d);
  theta0[3] ~ normal(f, g);
  
  
  Sigma ~ inv_wishart(4, D);
  
  for(i in 1:S){
    theta[i,] ~ multi_normal(theta0, Sigma);
  }
  
}

generated quantities{
  vector[3] SD;
  matrix[3, 3] CORR;
  for(i in 1:3){
    SD[i] = sqrt(Sigma[i, i]);
  }
  for(i in 1:3){
    for(j in 1:3){
      CORR[i, j] = Sigma[i, j] / (SD[i] * SD[j]);
    }
  }
  
}

