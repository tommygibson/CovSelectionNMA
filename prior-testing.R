###### Sampling from the prior of regularized horseshoe for covariance matrix
###### Regularized horseshoe only omega parameters
###### gamma parameters have normal priors

library(rstan)
library(tidyverse)
library(here)

options(mc.cores = parallel::detectCores())
###### Drawing from the prior


reghorse.prior.lambda.model <- stan_model(here("Models", "regularized-horseshoe-prior-lambda.stan"))

reghorse.prior.lambda.dat <- list(ntrt = 3, scale_gamma = 1 / sqrt(2),
                           scale_global_lambda = 1 / sqrt(2), nu_global_lambda = 3, nu_local_lambda = 1,
                           slab_scale_lambda = 4, slab_df_lambda = 1)

reghorse.prior.lambda.fit <- sampling(reghorse.prior.lambda.model, data = reghorse.prior.lambda.dat,
                               pars = c("SD", "CORR"), cores = 1)
reghorse.prior.lambda.sd <- as.data.frame(rstan::extract(reghorse.prior.lambda.fit, pars = "SD")$SD)
reghorse.prior.lambda.corr <- as.data.frame(rstan::extract(reghorse.prior.lambda.fit, pars = "CORR")$CORR)

names(reghorse.prior.lambda.sd) <- c("SD.1", "SD.2", "SD.3")
names(reghorse.prior.lambda.corr) <- c("CORR.1.1", "CORR.1.2", "CORR.1.3", 
                                "CORR.2.1", "CORR.2.2", "CORR.2.3", 
                                "CORR.3.1", "CORR.3.2", "CORR.3.3")

#PLOT PRIORS
par(mfrow = c(2, 2))
for(i in 1:3){
  hist(reghorse.prior.lambda.sd[reghorse.prior.lambda.sd[,i] < 5, i], breaks = 40)
}

par(mfrow = c(3, 3))
for(i in 1:ncol(reghorse.prior.lambda.corr)){
  hist(reghorse.prior.lambda.corr[, i], breaks = 40)
}

# see CDF values for sd < 1.5, < 0.1, < 0.01
apply(reghorse.prior.lambda.sd, 2, function(x){
  return(c(sum(x < 1.5) / length(x), sum(x < 0.1) / length(x), sum(x < 0.01) / length(x)))
})
apply(reghorse.prior.lambda.sd, 2, function(x){
  return(mean(x[x < 5]))
})

# see CDF values for |corr| < 0.5, < 0.1, < 0.01

apply(reghorse.prior.lambda.corr[,c(2, 3, 6)], 2, function(x){
  return(c(sum(abs(x) < 0.5) / length(x), sum(abs(x) < 0.1) / length(x), sum(abs(x) < 0.01) / length(x)))
})


################ 

reghorse.prior.model <- stan_model(here("Models", "regularized-horseshoe-prior.stan"))

reghorse.prior.dat <- list(ntrt = 3, 
                           scale_global_gamma = 2, scale_global_lambda = 1, 
                           nu_global_lambda = 1, nu_local_lambda = 1,
                           nu_global_gamma = 3, nu_local_gamma = 1,
                           slab_scale_lambda = 2, slab_df_lambda = 16,
                           slab_scale_gamma = 1, slab_df_gamma = 4)

reghorse.prior.fit <- sampling(reghorse.prior.model, data = reghorse.prior.dat,
                                      pars = c("SD", "CORR"), cores = 1)
reghorse.prior.sd <- as.data.frame(rstan::extract(reghorse.prior.fit, pars = "SD")$SD)
reghorse.prior.corr <- as.data.frame(rstan::extract(reghorse.prior.fit, pars = "CORR")$CORR)

names(reghorse.prior.sd) <- c("SD.1", "SD.2", "SD.3")
names(reghorse.prior.corr) <- c("CORR.1.1", "CORR.1.2", "CORR.1.3", 
                                       "CORR.2.1", "CORR.2.2", "CORR.2.3", 
                                       "CORR.3.1", "CORR.3.2", "CORR.3.3")

#PLOT PRIORS
par(mfrow = c(2, 2))
for(i in 1:3){
  hist(reghorse.prior.sd[reghorse.prior.sd[,i] < 5, i], breaks = 50)
}

par(mfrow = c(3, 3))
for(i in 1:ncol(reghorse.prior.corr)){
  hist(reghorse.prior.corr[, i], breaks = 50)
}

# see CDF values for sd < 1.5, < 0.1, < 0.01
apply(reghorse.prior.sd, 2, function(x){
  return(c(sum(x < 1.5) / length(x), sum(x < 0.1) / length(x), sum(x < 0.01) / length(x)))
})
apply(reghorse.prior.sd, 2, function(x){
  return(mean(x[x < 5]))
})

# see CDF values for |corr| < 0.5, < 0.1, < 0.01

apply(reghorse.prior.corr[,c(2, 3, 6)], 2, function(x){
  return(c(sum(abs(x) < 0.5) / length(x), sum(abs(x) < 0.1) / length(x), sum(abs(x) < 0.01) / length(x)))
})


#### For NMA

reghorse.prior.dat <- list(ntrt = 5, 
                           scale_global_gamma = 0.5, scale_global_lambda = 0.1, 
                           nu_global_lambda = 1, nu_local_lambda = 1,
                           nu_global_gamma = 3, nu_local_gamma = 1,
                           slab_scale_lambda = 4, slab_df_lambda = 4,
                           slab_scale_gamma = 4, slab_df_gamma = 4)

reghorse.prior.fit <- sampling(reghorse.prior.model, data = reghorse.prior.dat,
                               pars = c("SD", "CORR"), cores = 1)
reghorse.prior.sd <- as.data.frame(rstan::extract(reghorse.prior.fit, pars = "SD")$SD)
reghorse.prior.corr <- as.data.frame(rstan::extract(reghorse.prior.fit, pars = "CORR")$CORR)

#PLOT PRIORS
par(mfrow = c(2, 3))
for(i in 1:ncol(reghorse.prior.sd)){
  hist(reghorse.prior.sd[reghorse.prior.sd[,i] < 5, i], breaks = 40)
}

par(mfrow = c(5, 5))
for(i in 1:ncol(reghorse.prior.corr)){
  hist(reghorse.prior.corr[, i], breaks = 40)
}

# see CDF values for sd < 1.5, < 0.1, < 0.01
apply(reghorse.prior.sd, 2, function(x){
  return(c(sum(x < 1.5) / length(x), sum(x < 0.1) / length(x), sum(x < 0.01) / length(x)))
})
apply(reghorse.prior.sd, 2, function(x){
  return(mean(x[x < 5]))
})

# see CDF values for |corr| < 0.5, < 0.1, < 0.01

apply(reghorse.prior.corr[,c(2:5, 8:10, 14:15, 20)], 2, function(x){
  return(c(sum(abs(x) < 0.5) / length(x), sum(abs(x) < 0.1) / length(x), sum(abs(x) < 0.01) / length(x)))
})
