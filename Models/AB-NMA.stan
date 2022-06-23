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
  matrix[ntrt, ntrt] Lambda; // parameter lambda in Inverse - Wishart , generally indentity matrix
  real<lower = 1> nu; // parameter nu in Inverse - Wishart , generally ntrt +1
}
parameters {
  matrix [nstudy, ntrt] vi; // random effects
  vector [ntrt] mu; // fixed effect for treatment
  cov_matrix [ntrt] Sigma; // covariance matrix
}
transformed parameters {
  vector<lower = 0, upper = 1>[len] p; // Study_specified Absolute Risk
  for(i in 1:len){
    p[i] = inv_logit(mu[t[i]]+ vi[s[i], t[i]]);
  }
}
model {
  for(i in 1:len){
    r[i] ~ binomial (totaln[i], p[i]);
  }
  for(j in 1:nstudy){
    vi[j,]~ multi_normal(zeros, Sigma);
  }
  for(j in 1:ntrt){
    mu[j] ~ normal(0, 10);
  }
  Sigma ~ inv_wishart(nu, Lambda);
}
generated quantities {
  real AR[ntrt]; // Populaiton averaged absolute risk
  real MU[ntrt];
  real SD[ntrt]; // standard deviation of covariance matrix
  real OR[ntrt, ntrt]; // Odds Ratio
  real LOR[ntrt , ntrt]; // Log odds Ratio
  real cLOR [ntrt, ntrt]; // cLog odds Ratio
  real CORR [ntrt, ntrt]; // correlation matrix
  int rk[ntrt]; // rank of AR
  int rank_prob [ntrt, ntrt]; // rank probability

  real sumll = 0; // sum of log likelihood
  real log_likelihood [len];
  real likelihood [len];
  for(j in 1:ntrt){
    AR[j] = inv_logit(mu[j]/ sqrt(1 + Sigma [j,j] * 256 / 75 / pi() / pi()));
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
