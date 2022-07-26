# Testing stan version of 3RE vs jags

library(R2jags)
library(rstan)
library(tidyverse)
library(here)
library(MASS)

# generate data
options(mc.cores = parallel::detectCores())
set.seed(711)
S <- 10
N <- 300
beta0 <- nu0 <- log(.15 / .85)
delta0 <- 2
rho.beta.delta <- -0.7
rho.beta.nu <- 0
rho.delta.nu <- 0.7
sigma.beta <- 0.5
sigma.delta <- 1
sigma.nu <- 0.05

hyper.mean <- c(beta0, delta0, nu0)
hyper.var <- matrix(c(sigma.beta ^ 2, rho.beta.delta * sigma.beta * sigma.delta, rho.beta.nu * sigma.beta * sigma.nu,
                      rho.beta.delta * sigma.beta * sigma.delta, sigma.delta ^ 2, rho.delta.nu * sigma.delta * sigma.nu,
                      rho.beta.nu * sigma.beta * sigma.nu, rho.delta.nu * sigma.delta * sigma.nu, sigma.nu ^ 2), 
                    nrow = 3, byrow = TRUE)

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

re3.dat.jags <- list(S = S, y = y, n = n, n.tot = bigN,
                     a = -1, b = 0.25, c = 0, d = 0.25, e = -1, f = 0.25,
                     B.beta = 4, B.delta = 4, B.nu = 4)

re3.dat.stan <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                     a = -1, b = 2, c = 0, d = 2, f = -1, g = 2)
re3.dat.reghorse <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                         a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                         scale_global_lambda = 0.5, scale_global_gamma = 1,
                         nu_global_lambda = 1, nu_global_gamma = 1,
                         nu_local_lambda = 1, nu_local_gamma = 1,
                         slab_scale_lambda = 2, slab_df_lambda = 10,
                         slab_scale_gamma = 2, slab_df_gamma = 4)

re3.dat.lncass <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                       a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                       scale_local_lambda = 4, scale_local_gamma = 4,
                       tau_lambda = 1, tau_gamma = 1)

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
re3.lncass.model <- stan_model(here("Models", "3RE-LN-CASS.stan"))

# fit stan models
re3.stan <- sampling(re3.model, data = re3.dat.stan,
                     chains = 4, iter = 4000)#,
                     #control = list(adapt_delta = 0.95, max_treedepth = 15))
re3.lkj.stan <- sampling(re3.lkj.model, data = re3.dat.stan,
                         chains = 4, iter = 4000, #warmup = 4000,
                         pars = c("theta0", "Sigma"),
                         control = list(adapt_delta = .99))
re3.horseshoe.stan <- sampling(re3.horseshoe.model, data = re3.dat.stan,
                               chains = 4, iter = 8000,
                               pars = c("theta0", "Sigma"),
                               control = list(adapt_delta = 0.99, max_treedepth = 15),
                               cores = 4)
re3.reghorseshoe <- sampling(re3.reghorseshoe.model, data = re3.dat.reghorse,
                             chains = 4, iter = 4000, #warmup = 6000,
                             pars = c("theta0", "Sigma"),
                             control = list(adapt_delta = .99, max_treedepth = 15),
                             cores = 4)
re3.lncass <- sampling(re3.lncass.model, data = re3.dat.lncass,
                       chains = 4, iter = 4000,
                       pars = c("theta0", "Sigma"),
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       cores = 4)

round(re3.jags$BUGSoutput$summary, 4)
round(summary(re3.stan, pars = c("beta0", "delta0", "nu0", "sigma_beta", "sigma_delta", "sigma_nu"))$summary, 4)
round(summary(re3.lkj.stan)$summary, 4)
round(summary(re3.horseshoe.stan)$summary, 4)
round(summary(re3.reghorseshoe)$summary, 4)
round(summary(re3.lncass)$summary, 4)

