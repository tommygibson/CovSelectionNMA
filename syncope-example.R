##### 3RE model with covariance selection for syncope data

library(rstan)
library(tidyverse)
library(here)
library(loo)

source(here("CTS-functions.R"))


# define the models

re3.std.model <- stan_model(here("Models", "3RE.stan"))
re3.iw.model <- stan_model(here("Models", "3RE-IW.stan"))
re3.lkj.model <- stan_model(here("Models", "3RE-LKJ.stan"))
re3.reghorse.model <- stan_model(here("Models", "3RE-regularized-horseshoe-cholesky.stan"))
re3.reglambda.model <- stan_model(here("Models", "3RE-regularized-lambda-cholesky.stan"))

# load syncope data, filter out vars with only 2-4 observations
# also remove studies without counts

syncope <- read_csv(here("syncope.cleaned.csv")) %>%
  group_by(Variable) %>%
  filter(!is.na(n_i0),
         n() > 4)

# start with male gender as an example
male <- syncope %>%
  ungroup %>%
  filter(Variable == "Male Gender")

y <- as.matrix(male %>% 
  dplyr::select(y_i1, y_i0) )
n <- as.matrix(male %>% 
  dplyr::select(n_i1, n_i0))
bigN <- rowSums(n)
S <- nrow(y)

nu.hat <- log((n[,1] / bigN) / (1 - n[,1] / bigN))
beta.hat <- apply(log((y/n) / (1 - y/n)), 1, mean)
delta.hat <- male$lnORhat

sd.nu.hat <- sd(nu.hat)
sd.beta.hat <- sd(beta.hat)
sd.delta.hat <- sd(delta.hat)

sd.sd.nu.hat <- sd.nu.hat * sqrt(1 - 2 / (S - 1) * (gamma(S / 2) / gamma((S - 1) / 2))^2)
sd.sd.beta.hat <- sd.beta.hat * sqrt(1 - 2 / (S - 1) * (gamma(S / 2) / gamma((S - 1) / 2))^2)
sd.sd.delta.hat <- sd.delta.hat * sqrt(1 - 2 / (S - 1) * (gamma(S / 2) / gamma((S - 1) / 2))^2)

(corhat.1.2 <- cor(beta.hat, delta.hat))
(corhat.1.3 <- cor(beta.hat, nu.hat))
(corhat.2.3 <- cor(delta.hat, nu.hat))

mean.sd.sd <- mean(c(sd.sd.nu.hat, sd.sd.beta.hat, sd.sd.delta.hat))
scale_global_lambda <- 2 / (3 - 2) * mean.sd.sd

re3.dat.stan <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                     a = 0, b = 2, c = 0, d = 2, f = 0, g = 2)
re3.dat.reghorse <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                         a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                         scale_global_lambda = scale_global_lambda, scale_global_gamma = 2,
                         nu_global_lambda = 1, nu_global_gamma = 3,
                         nu_local_lambda = 1, nu_local_gamma = 1,
                         slab_scale_lambda = 2, slab_df_lambda = 4,
                         slab_scale_gamma = 2, slab_df_gamma = 4)
re3.dat.reglambda <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                         a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                         scale_global_lambda = scale_global_lambda, scale_gamma = 8,
                         nu_global_lambda = 1,  nu_local_lambda = 1,
                         slab_scale_lambda = 2, slab_df_lambda = 4)
# scale_global_lambda = 1, scale_global_gamma = 1,
# nu_global_lambda = 1, nu_global_gamma = 1,
# nu_local_lambda = 1, nu_local_gamma = 1,
# slab_scale_lambda = 2, slab_df_lambda = 4,
# slab_scale_gamma = 1, slab_df_gamma = 8
re3.stan <- sampling(re3.model, data = re3.dat.stan,
                     chains = 4, iter = 4000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15))
re3.iw.stan <- sampling(re3.iw.model, data = re3.dat.stan,
                        chains = 4, iter = 4000,
                        pars = c("theta0", "SD", "CORR"))#,
                        #control = list(adapt_delta = 0.95))

re3.lkj.stan <- sampling(re3.lkj.model, data = re3.dat.stan,
                         chains = 4, iter = 4000, # iter = 10000, warmup = 6000,
                         pars = c("theta0", "Sigma"),
                         control = list(adapt_delta = .99, max_treedepth = 15),
                         cores = 4)

re3.reghorse <- sampling(re3.reghorsechol.model, data = re3.dat.reghorse,
                             chains = 4, iter = 5000, warmup = 2500,
                             pars = c("theta0", "SD", "CORR"),
                             control = list(adapt_delta = .99, max_treedepth = 20),
                             cores = 4)
re3.reglambda <- sampling(re3.reglambda.model, data = re3.dat.reglambda,
                          chains = 4, iter = 12000, warmup = 6000,
                          pars = c("theta0", "Sigma"),
                          control = list(adapt_delta = 0.99, max_treedepth = 15),
                          cores = 4)

round(summary(re3.stan, pars = c("beta0", "delta0", "nu0", "sigma_beta", "sigma_delta", "sigma_nu", "lp__"))$summary, 4)
round(summary(re3.iw.stan)$summary, 4)
round(summary(re3.lkj.stan)$summary, 4)
round(summary(re3.reghorse)$summary, 4)
round(summary(re3.reglambda)$summary, 4)

# reglambda.sims <- rstan::extract(re3.reglambda, pars = c("theta0", "Sigma"))
make_mean_var_sims <- function(sims.list){
  theta <- sims.list$theta0
  Sigma.list <- list()
  for(i in 1:dim(sims.list$Sigma)[1]){
    Sigma.list[[i]] <- as.vector(sims.list$Sigma[i,,])
  }
  Sigma <- do.call(rbind, Sigma.list)
  return(cbind(theta, Sigma))
}

reglambda.sims <- make_mean_var_sims(rstan::extract(re3.reglambda, pars = c("theta0", "Sigma")))
lkj.sims <- make_mean_var_sims(rstan::extract(re3.lkj.stan, pars = c("theta0", "Sigma")))
cts0_from_one_posterior(reglambda.sims[1,])
reglambda.cts0.summary <- as.data.frame(cts_summary_from_gamma(reglambda.sims, M = 1000)) %>%
  type_convert() %>%
  rename(CTS =  V1,
         Mean = V2,
         SD =   V3)
lkj.cts0.summary <- as.data.frame(cts_summary_from_gamma(lkj.sims, M = 1000)) %>%
  type_convert() %>%
  rename(CTS = V1,
         Mean = V2,
         SD = V3)
