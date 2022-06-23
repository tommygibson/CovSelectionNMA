### testing AB-NMA models on sample datasets

library(rstan)
library(R2jags)
library(metadat)
library(tidyverse)
library(here)
library(MASS)


rstan_options("auto_write" = TRUE)
options(mc.cores = 4)


#### Generating a test dataset
set.seed(616)
N <- 300
ntrt <- 5
nstudy <- 40

# absolute mean effects for each treatment arm
mu <- c(-1.75, -2, -2.5, -1.5, -2.75)
# full mean matrix, because theoretically each study contained each arm and they're MCAR
mu.tot <- matrix(rep(mu, nstudy), nrow = nstudy, byrow = TRUE)

# construct covariance matrix of random effects with the Lambda %*% Gamma %*% t(Gamma) %*% Lambda formulation
Lambda <- diag(c(0.3, 0.1, 0.4, 0.5, 0.2))
gamma <- c(0.1, 0.4, 0.6, 0, 0, 0, 0.8, 0.2, 0.1, 0.2)
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
v.tot <- mvrnorm(nstudy, rep(0, ntrt), Sigma = Sigma)

p.tot <- 1 / (1 + exp(-(mu.tot + v.tot)))

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
for(i in 2:ntrt){
  non.missing.pattern[, i] <- c(rep(0, 10 * (i - 2)), rep(1, 10), rep(0, 10 * (5 - i)))
}


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

test.dat.stan <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                      s = s, r = r, totaln = totaln,
                      higher_better = 0, zeros = zeros,
                      Lambda = Lambda, 
                      nu = nu)
dat.stan.hhc <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                     s = s, r = r, totaln = totaln,
                     higher_better = 0, zeros = zeros,
                     Lambda = Lambda, 
                     nu = nu, uni_upper = 3)
dat.stan.lkj <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                     s = s, r = r, totaln = totaln,
                     higher_better = 0, zeros = zeros,
                     Lambda = Lambda, 
                     nu = nu, uni_upper = 3, eta = 1.5)

jags.dat.select <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                        s = s, r = r, totaln = totaln, 
                        zeros = zeros, c = c, pi = pi)

ab.nma.select.inits <- function(){
  list(mu = runif(ntrt, -1, 1),
       vi = matrix(runif(ntrt * nstudy, -1, 1), nrow = nstudy),
       gamma.z = rbinom((ntrt * (ntrt - 1) / 2), 1, .5),
       gamma.n = runif((ntrt * (ntrt - 1) / 2), -1, 1),
       lambda = runif(ntrt, 0.1, 0.5))
}

ab.nma.stan <- stan_model(here("Models", "AB-NMA.stan"))
stan.hhc <- stan_model(here("Models", "AB-NMA-HHC.stan"))
stan.lkj <- stan_model(here("Models", "AB-NMA-HHC-LKJ.stan"))

fit.stan.IW <- sampling(ab.nma.stan, pars = c("AR", "CORR", "LOR", "MU"), data = test.dat.stan, 
                       iter = 4000, cores = 4, thin = 4)
fit.stan.HHC <- sampling(stan.hhc, pars = c("AR", "MU"), data = dat.stan.hhc,
                         iter = 4000, cores = 4, thin = 4,
                         control = list(adapt_delta = 0.99, max_treedepth = 15))
fit.stan.LKJ <- sampling(stan.lkj, pars = c("AR", "MU"), data = dat.stan.lkj,
                         iter = 4000, cores = 4, thin = 4,
                         control = list(adapt_delta = 0.99, max_treedepth = 15))
ab.nma.params <- c("AR", "CORR", "LOR", "MU")

ab.select.fit <- do.call(jags.parallel,
                         list(names(jags.dat.select), parameters.to.save = ab.nma.params,
                              inits = ab.nma.select.inits,
                              model.file = here("Models", "AB-NMA-select.txt"), DIC = FALSE,
                              n.iter = 40000, n.burnin = 10000, n.chains = 4, n.thin = 12))



########
nma.dat.1 <- dat.dong2013
len <- dim(nma.dat.1)[1]
ntrt <- length(unique(nma.dat.1$treatment))
nstudy <- max(nma.dat.1$id)
t <- as.numeric(as.factor(nma.dat.1$treatment))
s <- nma.dat.1$id
r <- nma.dat.1$death
totaln <- nma.dat.1$randomized
zeros <- rep(0, ntrt)
Lambda <- diag(ntrt)
nu <- ntrt + 1


stan.dat.1 <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                   s = s, r = r, totaln = totaln,
                   higher_better = 0, zeros = zeros,
                   Lambda = Lambda, 
                   nu = nu)
jags.dat.1 <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                   s = s, r = r, totaln = totaln, 
                   zeros = zeros, Lambda = Lambda, nu = nu, pi = pi)

ab.nma.inits <- function(){
  list(mu = runif(ntrt, -1, 1),
       vi = matrix(runif(ntrt * nstudy, -1, 1), nrow = nstudy))
}

ab.nma.stan <- stan_model(here("Models", "AB-NMA.stan"))
params <- c("LOR", "CORR")


fit.1.stan <- sampling(ab.nma.stan, pars = c("LOR", "CORR"), data = stan.dat.1, 
                  iter = 8000, cores = 4, thin = 4)
fit.1.jags <- do.call(jags.parallel,
                      list(names(jags.dat.1), ab.nma.inits, params, 
                   model.file = here("Models", "AB-NMA.txt"),
                   n.iter = 80000, n.burnin = 10000, n.chains = 4, n.thin = 8, DIC = FALSE))

summary(fit.1.stan, pars = c("Sigma"))$summary
fit.1.jags$BUGSoutput$summary


## testing ifelse in jags
ntrt <- 4
ifelse.dat <- list(ntrt = ntrt, p = 0.5)
ifelse.params <- c("Sigma")
ifelse.inits <- function(){
  list(
    gamma.z = rbinom(6, 1, 0.5),
    gamma.n = rnorm(6),
    lambda = rnorm(ntrt)
  )
}

ifelse.fit <- jags(ifelse.dat, ifelse.inits, ifelse.params,
                   model.file = here("Models", "ifelse.test.txt"), DIC = FALSE)

# fitting AB-NMA-select to simple dataset...
len <- dim(nma.dat.1)[1]
ntrt <- length(unique(nma.dat.1$treatment))
nstudy <- max(nma.dat.1$id)
t <- as.numeric(as.factor(nma.dat.1$treatment))
s <- nma.dat.1$id
r <- nma.dat.1$death
totaln <- nma.dat.1$randomized
zeros <- rep(0, ntrt)
c <- 0.5

jags.dat.select <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                   s = s, r = r, totaln = totaln, 
                   zeros = zeros, c = c, pi = pi)

ab.nma.select.inits <- function(){
  list(mu = runif(ntrt, -1, 1),
       vi = matrix(runif(ntrt * nstudy, -1, 1), nrow = nstudy),
       gamma.z = rbinom((ntrt * (ntrt - 1) / 2), 1, .5),
       gamma.n = runif((ntrt * (ntrt - 1) / 2), -1, 1),
       lambda = runif(ntrt, 0.1, 0.5))
}
ab.nma.params <- c("AR", "CORR", "LOR", "MU")

ab.select.fit <- do.call(jags.parallel,
                         list(names(jags.dat.select), parameters.to.save = ab.nma.params,
                         inits = ab.nma.select.inits,
                         model.file = here("Models", "AB-NMA-select.txt"), DIC = FALSE,
                         n.iter = 250000, n.burnin = 200000, n.chains = 4, n.thin = 20))

# compare stan fit and selection model fit
# diffs in sds first
LOR.mean.stan <- matrix(round(summary(fit.1.stan, pars = c("LOR"))$summary[,1], 3), nrow = 6, byrow =F)
LOR.mean.jags <- matrix(round(ab.select.fit$BUGSoutput$summary[37:72,2], 3), nrow = 6, byrow = T)
LOR.sd.stan <- matrix(round(summary(fit.1.stan, pars = c("LOR"))$summary[,3], 3), nrow = 6, byrow = F)
LOR.sd.select <- matrix(round(ab.select.fit$BUGSoutput$summary[37:72,2], 3), nrow = 6, byrow = T)

LOR.sd.stan - LOR.sd.select
