// 3RE model WITHOUT covariance selection
// standard model to be compared with JAGS version



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
  vector[S] beta;
  vector[S] delta;
  vector[S] nu;
  real beta0;
  real delta0;
  real nu0;
  real<lower = 0> sigma_beta;
  real<lower = 0> sigma_delta;
  real<lower = 0> sigma_nu;
}

transformed parameters {
  vector[S] pi1;
  vector[S] pi0;
  vector[S] psi;
  for(i in 1:S){
    pi1[i] = inv_logit(beta[i] + delta[i] / 2);
    pi0[i] = inv_logit(beta[i] - delta[i] / 2);
    psi[i] = inv_logit(nu[i]);
  }
}

model {
  // binomial likelihood
  y1 ~ binomial(n1, pi1);
  y0 ~ binomial(n0, pi0);
  n1 ~ binomial(N, psi);
  
  // priors
  beta ~ normal(beta0, sigma_beta);
  delta ~ normal(delta0, sigma_delta);
  nu ~ normal(nu0, sigma_nu);
  
  // hyperpriors
  beta0 ~ normal(a, b);
  delta0 ~ normal(c, d);
  nu0 ~ normal(f, g);
  
  sigma_beta ~ cauchy(0, 0.75);
  sigma_delta ~ cauchy(0, 0.75);
  sigma_nu ~ cauchy(0, 0.75);
  
}

