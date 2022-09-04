####### Simulation to test the added study-specific random effect AB-NMA model
library(rstan)
library(tidyverse)
library(here)
library(MASS)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
set.seed(825)

nma.std <- stan_model(here("Models", "AB-NMA.stan"))
#nma.hhc <- stan_model(here("Models", "AB-NMA-HHC.stan"))
nma.lkj <- stan_model(here("Models", "AB-NMA-LKJ.stan"))
#nma.reghorse <- stan_model(here("Models", "AB-NMA-regularized-horseshoe.stan"))
nma.reglambda <- stan_model(here("Models", "AB-NMA-regularized-lambda.stan"))
nma.reglambda.sepvar <- stan_model(here("Models", "AB-NMA-sep-var-reglambda.stan"))
nma.reglambda.sepvar.hier <- stan_model(here("Models", "AB-NMA-sepvar-hier.stan"))
nma.reglambda.sepvar.hier.freegamma <- stan_model(here("Models", "AB-NMA-sepvar-hier-freegamma.stan"))

nstudy <- 40
N <- 400
ntrt <- 5

mu <- c(-2, -2.5, -2.25, -2.75, -1.75)
mu.tot <- matrix(rep(mu, nstudy), nrow = nstudy, byrow = TRUE)

beta <- rnorm(S, 0, 0.5)

# construct covariance matrix of random effects with the Lambda %*% Gamma %*% t(Gamma) %*% Lambda formulation
Lambda <- diag(c(0.3, 0.001, 0.8, 0.5, 0.2))#, .001))
gamma <- c(0.6, 
           -0.2, 0.8, 
           0, 0, 0,
           0.8, -0.2, 0, 0.2,
           0, 0, 0, 0, 0)
Gamma <- diag(ntrt)
index <- 1
for(i in 2:ntrt){
  Gamma[i, 1:(i - 1)] <- gamma[index:(index + i - 2)]
  index <- index + (i - 1)
}

Sigma <- Lambda %*% Gamma %*% t(Gamma) %*% Lambda

# let's see the correlation mat just for kicks
P <- diag(ntrt)
for(i in 1:ntrt){
  for(j in 1:(ntrt - 1)){
    P[i, j] <- ifelse(i == j, 1, Sigma[i, j] / sqrt(Sigma[i, i] * Sigma[j, j]))
    P[j, i] <- P[i, j]
  }
}

# random effects, multivariate normal
v.tot <- matrix(nrow = nstudy, ncol = ntrt)
for(i in 1:nstudy){
  v.tot[i,] <- mvrnorm(1, mu + beta[i], Sigma = Sigma)
}


p.tot <- 1 / (1 + exp(-(v.tot)))

ri <- matrix(NA, nrow = nstudy, ncol = ntrt)
for(i in 1:nstudy){
  for(j in 1:ntrt){
    ri[i, j] <- rbinom(1, size = N, prob = p.tot[i, j])
  }
}

# decide which studies have missing values
# for now we'll say that there are only 2-arm studies
# each arm (other than "placebo" (1) gets 40 / 4 = 10 studies non-missing)
# we'll come up with something better later
non.missing.pattern <- matrix(0, nrow = nstudy, ncol = ntrt)
non.missing.pattern[,1] <- 1
non.missing.pattern[,2] <- c(rep(1, 12), rep(0, 28))
non.missing.pattern[,3] <- c(rep(0, 10), rep(1, 12), rep(0, 18))
non.missing.pattern[,4] <- c(rep(0, 20), rep(1, 12), rep(0, 8))
non.missing.pattern[,5] <- c(rep(0, 28), rep(1, 12))


# alternative approach, randomize the missingness
len <- sum(non.missing.pattern)
t <- rep(0, len)
s <- rep(0, len)
r <- rep(0, len)
totaln <- rep(N, len)
index <- 1
for(i in 1:nstudy){
  for(j in 1:ntrt){
    if(non.missing.pattern[i, j] == 1){
      t[index] <- j
      s[index] <- i
      r[index] <- ri[i, j]
      index <- index + 1
    }
  }
}

zeros <- rep(0, ntrt)
Lambda <- diag(ntrt)
nu <- ntrt + 1
c <- 0.5

std.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                      s = s, r = r, totaln = totaln,
                      higher_better = 0, zeros = zeros,
                      Lambda = Lambda, 
                      nu = nu)
dat.stan.hhc <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                     s = s, r = r, totaln = totaln,
                     higher_better = 0, zeros = zeros,
                     Lambda = Lambda, 
                     nu = nu, uni_upper = 3)
lkj.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                     s = s, r = r, totaln = totaln,
                     higher_better = 0, zeros = zeros,
                     Lambda = Lambda, 
                     nu = nu, eta = 1)
reghorse.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                     s = s, r = r, totaln = totaln,
                     higher_better = 0, zeros = zeros,
                     scale_global_lambda = 0.1, scale_global_gamma = 0.4,
                     nu_global_lambda = 1, nu_global_gamma = 1,
                     nu_local_lambda = 1, nu_local_gamma = 1,
                     slab_scale_lambda = 4, slab_scale_gamma = 2,
                     slab_df_lambda = 1, slab_df_gamma = 1)

reglambda.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                      s = s, r = r, totaln = totaln,
                      higher_better = 0, zeros = zeros,
                      scale_gamma = 2,
                      scale_global_lambda = 0.1, nu_global_lambda = 3, nu_local_lambda = 1,
                      slab_scale_lambda = 2, slab_df_lambda = 1)

reglambda.sepvar.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                             s = s, r = r, totaln = totaln,
                             higher_better = 0, zeros = zeros,
                             scale_gamma = 4,
                             scale_global_lambda = 0.1, nu_global_lambda = 1, nu_local_lambda = 1,
                             slab_scale_lambda = 2, slab_df_lambda = 4)

reglambda.sepvar.freegamma.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                                       s = s, r = r, totaln = totaln,
                                       higher_better = 0, zeros = zeros,
                                       scale_global_lambda = 0.1, nu_global_lambda = 1, nu_local_lambda = 1,
                                       slab_scale_lambda = 2, slab_df_lambda = 4)



fit.1.stan <- sampling(nma.std, pars = c("SD", "CORR", "MU", "AR",
                                         "LOR", "log_likelihood"), 
                       data = std.dat, iter = 1000, cores = 4,
                       control = list(adapt_delta = 0.99))

fit.1.lkj <- sampling(nma.lkj, pars = c("SD", "CORR", "MU", "AR",
                                        "LOR", "log_likelihood"), 
                      data = lkj.dat, iter = 1000, cores = 4,
                      control = list(adapt_delta = 0.99, max_treedepth = 15))

fit.1.reglambda <- sampling(nma.reglambda, 
                            pars = c("SD", "CORR", "MU", "AR",
                                     "LOR", "log_likelihood"), 
                            data = reglambda.dat, iter = 1000, #warmup = 4000, cores = 4,
                            control = list(adapt_delta = 0.99, max_treedepth = 15))
# fit.1.reglambda.sepvar <- sampling(nma.reglambda.sepvar,
#                                    pars = c("SD", "sigma_beta", "CORR",
#                                             "MU", "AR", "LOR", "log_likelihood"),
#                                    data = reglambda.dat,
#                                    iter = 12000, cores = 4,
#                                    control = list(adapt_delta = 0.99, max_treedepth = 15))
fit.1.sepvar.hier <- sampling(nma.reglambda.sepvar.hier,
                              pars = c("SD", "sigma_beta", "CORR",
                                       "MU", "AR", "LOR", "log_likelihood"),
                              data = reglambda.sepvar.dat,
                              iter = 1000, cores = 4,
                              control = list(adapt_delta = 0.99, max_treedepth = 15))
fit.1.sepvar.freegamma <- sampling(nma.reglambda.sepvar.hier.freegamma,
                                   pars = c("SD", "sigma_beta", "CORR",
                                            "MU", "AR", "LOR", "log_likelihood"),
                                   data = reglambda.sepvar.freegamma.dat,
                                   iter = 1000, cores = 4,
                                   control = list(adapt_delta = 0.99, max_treedepth = 15))
