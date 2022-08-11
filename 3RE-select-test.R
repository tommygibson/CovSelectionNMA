# Testing stan version of 3RE vs jags

library(R2jags)
library(rstan)
library(tidyverse)
library(here)
library(MASS)

options(mc.cores = parallel::detectCores())

# generate data
set.seed(802)
S <- 7
N <- 400
beta0 <- nu0 <- log(.15 / .85)
delta0 <- 2
rho.beta.delta <- 0.7
rho.beta.nu <- 0
rho.delta.nu <- 0.7
sigma.beta <- 0.5
sigma.delta <- 1
sigma.nu <- 0.01

Lambda <- diag(c(.5, .5, 0.025))
Gamma <- diag(3)
Gamma[2,1] <- 0.6
Gamma[3,1] <- 0
Gamma[3,2] <- -0.6
hyper.var <- Lambda %*% Gamma %*% t(Gamma) %*% Lambda
cormat <- matrix(nrow = 3, ncol = 3)
for(i in 1:3){
  for(j in 1:3){
    cormat[i, j] <- hyper.var[i, j] / (sqrt(hyper.var[i, i]) * sqrt(hyper.var[j, j]))
  }
}

hyper.mean <- c(beta0, delta0, nu0)
# hyper.var <- matrix(c(sigma.beta ^ 2, rho.beta.delta * sigma.beta * sigma.delta, rho.beta.nu * sigma.beta * sigma.nu,
#                       rho.beta.delta * sigma.beta * sigma.delta, sigma.delta ^ 2, rho.delta.nu * sigma.delta * sigma.nu,
#                       rho.beta.nu * sigma.beta * sigma.nu, rho.delta.nu * sigma.delta * sigma.nu, sigma.nu ^ 2), 
#                     nrow = 3, byrow = TRUE)

theta <- mvrnorm(S, hyper.mean, hyper.var)
beta <- theta[,1]
delta <- theta[,2]
nu <- theta[,3]

pi1 <- 1 / (1 + exp(-(beta + delta / 2)))
pi0 <- 1 / (1 + exp(-(beta - delta / 2)))
psi <- 1 / (1 + exp(-nu))

pi11 <- pi1 * psi
pi10 <- (1 - pi1) * psi
pi01 <- pi0 * (1 - psi)
pi00 <- (1 - pi0) * (1 - psi)

probs <- cbind(pi11, pi10, pi01, pi00)

tabs <- t(apply(probs, 1, function(x) {
  rmultinom(1, N, prob = x)
}))

y <- tabs[,c(1, 3)]
n <- cbind(rowSums(tabs[, c(1, 2)]),
           rowSums(tabs[, c(3, 4)]))
bigN <- rep(N, S)

nu.hat <- log((n[,1] / bigN) / (1 - n[,1] / bigN))
beta.hat <- apply(log((y/n) / (1 - y/n)), 1, mean)
delta.hat <- log((tabs[,1] * tabs[,4]) / (tabs[,2] * tabs[,3]))

sd.nu.hat <- sd(nu.hat)
sd.beta.hat <- sd(beta.hat)
sd.delta.hat <- sd(delta.hat)

sd.sd.nu.hat <- sd.nu.hat * sqrt(1 - 2 / (S - 1) * (gamma(S / 2) / gamma((S - 1) / 2))^2)
sd.sd.beta.hat <- sd.beta.hat * sqrt(1 - 2 / (S - 1) * (gamma(S / 2) / gamma((S - 1) / 2))^2)
sd.sd.delta.hat <- sd.delta.hat * sqrt(1 - 2 / (S - 1) * (gamma(S / 2) / gamma((S - 1) / 2))^2)

(corhat.1.2 <- cor(beta.hat, delta.hat))
(corhat.1.3 <- cor(beta.hat, nu.hat))
(corhat.2.3 <- cor(delta.hat, nu.hat))

re3.dat.jags <- list(S = S, y = y, n = n, n.tot = bigN,
                     a = -1, b = 0.25, c = 0, d = 0.25, e = -1, f = 0.25,
                     B.beta = 4, B.delta = 4, B.nu = 4)

re3.dat.stan <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                     a = 0, b = 2, c = 0, d = 2, f = 0, g = 2)
re3.dat.reghorse <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                         a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                         scale_global_lambda = 0.5, scale_global_gamma = 1,
                         nu_global_lambda = 1, nu_global_gamma = 1,
                         nu_local_lambda = 1, nu_local_gamma = 1,
                         slab_scale_lambda = 1, slab_df_lambda = 1,
                         slab_scale_gamma = 0.5, slab_df_gamma = 1)
re3.dat.reglambda <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                          a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                          scale_gamma = 1, scale_global_lambda = 0.5, 
                          nu_global_lambda = 1, nu_local_lambda = 1,
                          slab_scale_lambda = 1, slab_df_lambda = 4)
# 
# re3.dat.lncass <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
#                        a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
#                        scale_local_lambda = 10, scale_local_gamma = 10,
#                        tau_lambda = 0.5, tau_gamma = 0.5)

init.gen.3re <- function(){
  list(beta0 = runif(1, -1, 1),
       delta0 = runif(1, -1, 1),
       nu0 = runif(1, -1, 1),
       beta = runif(S, -1, 1),
       delta = runif(S, -1, 1),
       nu = runif(S, -1, 1),
       sigma.beta = runif(1, 0.2, 1),
       sigma.delta = runif(1, 0.2, 1),
       sigma.nu = runif(1, 0.2, 1)
       )
}

re3.params.jags <- c("beta0", "delta0", "nu0", "sigma.beta", "sigma.delta", "sigma.nu")

re3.jags <- jags(re3.dat.jags, init.gen.3re, re3.params.jags, 
                 here("Models", "3RE.txt"),
                 n.chains = 4, n.iter = 10000, n.burnin = 4000, n.thin = 2, DIC = FALSE)
# compile stan models
re3.model <- stan_model(here("Models", "3RE.stan"))
re3.lkj.model <- stan_model(here("Models", "3RE-LKJ.stan"))
re3.horseshoe.model <- stan_model(here("Models", "3RE-horseshoe.stan"))
re3.reghorseshoe.model <- stan_model(here("Models", "3RE-regularized-horseshoe.stan"))
re3.reghorsechol.model <- stan_model(here("Models", "3RE-regularized-horseshoe-cholesky.stan"))
re3.reglambdachol.model <- stan_model(here("Models", "3RE-regularized-lambda-cholesky.stan"))
# ln-cass model is basically a failure
# re3.lncass.model <- stan_model(here("Models", "3RE-LN-CASS.stan"))

# fit stan models
re3.stan <- sampling(re3.model, data = re3.dat.stan,
                     chains = 4, iter = 4000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15))
re3.lkj.stan <- sampling(re3.lkj.model, data = re3.dat.stan,
                         chains = 4, iter = 4000, #warmup = 4000,
                         pars = c("theta0", "SD", "CORR"),
                         control = list(adapt_delta = .99, max_treedepth = 15))
re3.horseshoe.stan <- sampling(re3.horseshoe.model, data = re3.dat.stan,
                               chains = 4, iter = 8000,
                               pars = c("theta0", "Sigma", "corr"),
                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                               cores = 4)
re3.reghorseshoe <- sampling(re3.reghorseshoe.model, data = re3.dat.reghorse,
                             chains = 4, iter = 4000, #warmup = 6000,
                             pars = c("theta0", "SD", "CORR"),
                             control = list(adapt_delta = .99, max_treedepth = 15),
                             cores = 4)
re3.reghorsechol <- sampling(re3.reghorsechol.model, data = re3.dat.reghorse,
                             chains = 4, iter = 4000, #warmup = 6000,
                             pars = c("theta0", "SD", "CORR"),
                             control = list(adapt_delta = .99, max_treedepth = 15),
                             cores = 4)
re3.reglambdachol <- sampling(re3.reglambdachol.model, data = re3.dat.reglambda,
                              chains = 4, iter = 4000, #warmup = 6000,
                              pars = c("theta0", "SD", "CORR"),
                              control = list(adapt_delta = 0.99, max_treedepth = 15),
                              cores = 4)
# re3.lncass <- sampling(re3.lncass.model, data = re3.dat.lncass,
#                        chains = 4, iter = 4000,
#                        pars = c("theta0", "Sigma"),
#                        control = list(adapt_delta = 0.99, max_treedepth = 15),
#                        cores = 4)

round(re3.jags$BUGSoutput$summary, 4)
round(summary(re3.stan, pars = c("beta0", "delta0", "nu0", "sigma_beta", "sigma_delta", "sigma_nu", "lp__"))$summary, 4)
round(summary(re3.lkj.stan)$summary, 4)
round(summary(re3.horseshoe.stan)$summary, 4)
round(summary(re3.reghorseshoe)$summary, 4)
round(summary(re3.reghorsechol)$summary, 4)
round(summary(re3.reglambdachol)$summary, 4)



########## Crafting an example by hand with zero variance for each RE

set.seed(803)
S <- 20
N <- round(seq(300, 900, length.out = S))

n <- round(N %*% t(c(3/20, 17/20)))
bigN <- rowSums(n)

pi0 <- 0.07
pi1 <- runif(S, .2, .4)

y <- matrix(nrow = S, ncol = 2)
y[,1] <- rbinom(S, n[,1], pi1)
y[,2] <- rbinom(S, n[,2], pi0)

nu.hat <- log((n[,1] / bigN) / (1 - n[,1] / bigN))
beta.hat <- apply(log((y/n) / (1 - y/n)), 1, mean)
delta.hat <- log((y[,1] * (n[,2] - y[,2])) / ((n[,1] - y[,1]) * y[,2]))

re3.dat.stan <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                     a = 0, b = 2, c = 0, d = 2, f = 0, g = 2)
re3.dat.reghorse <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                         a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                         scale_global_lambda = 1, scale_global_gamma = 2,
                         nu_global_lambda = 1, nu_global_gamma = 3,
                         nu_local_lambda = 1, nu_local_gamma = 1,
                         slab_scale_lambda = 2, slab_df_lambda = 1,
                         slab_scale_gamma = 1, slab_df_gamma = 4)
re3.dat.reglambda <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                          a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                          scale_gamma = 1 / sqrt(2),
                          scale_global_lambda = 1, nu_global_lambda = 3, nu_local_lambda = 1,
                          slab_scale_lambda = 2, slab_df_lambda = 1)

re3.stan <- sampling(re3.model, data = re3.dat.stan,
                     chains = 4, iter = 4000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15))
# re3.iw <- samplign(re3.)
re3.lkj.stan <- sampling(re3.lkj.model, data = re3.dat.stan,
                         chains = 4, iter = 4000, # iter = 10000, warmup = 6000,
                         pars = c("theta0", "SD", "CORR"),
                         control = list(adapt_delta = .99, max_treedepth = 15),
                         cores = 4)
re3.horse <- sampling(re3.horseshoe.model, data = re3.dat.stan,
                      chains = 4, iter = 4000, 
                      control = list(adapt_delta = 0.99, max_treedepth = 15),
                      cores = 4)
# re3.reghorseshoe <- sampling(re3.reghorseshoe.model, data = re3.dat.reghorse,
#                              chains = 4, iter = 4000, #warmup = 6000,
#                              pars = c("theta0", "SD", "CORR"),
#                              control = list(adapt_delta = .99, max_treedepth = 15),
#                              cores = 4)
re3.reghorsechol <- sampling(re3.reghorsechol.model, data = re3.dat.reghorse,
                             chains = 4, iter = 8000, warmup = 4000,
                             pars = c("theta0", "SD", "CORR"),
                             control = list(adapt_delta = .99, max_treedepth = 20),
                             cores = 4)
# re3.reglambdachol <- sampling(re3.reglambdachol.model, data = re3.dat.reglambda,
#                               chains = 4, iter = 10000, warmup = 6000,
#                               pars = c("theta0", "SD", "CORR"),
#                               control = list(max_treedepth = 15),
#                               cores = 4)

## function to plot posterior from 3RE

lkj.samp <- as.data.frame(rstan::extract(re3.lkj.stan, pars = c("theta0", "SD", "CORR")))[, c(1:6, 8:9, 12)]
horse.samp <- as.data.frame(rstan::extract(re3.reghorsechol, pars = c("theta0", "SD", "CORR")))[, c(1:6, 8:9, 12)]
par(mfrow = c(3, 3))
for(i in 1:ncol(lkj.samp)){
  hist(lkj.samp[,i], breaks = 30)
}
for(i in 1:ncol(horse.samp)){
  hist(horse.samp[,i], breaks = 30)
}


round(summary(re3.stan, pars = c("beta0", "delta0", "nu0", "sigma_beta", "sigma_delta", "sigma_nu", "lp__"))$summary, 4)
round(summary(re3.lkj.stan)$summary, 4)
round(summary(re3.reghorseshoe)$summary, 4)
round(summary(re3.reghorsechol)$summary, 4)
round(summary(re3.reglambdachol)$summary, 4)
