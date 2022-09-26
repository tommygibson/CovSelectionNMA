### testing AB-NMA models on sample datasets

library(rstan)
library(tidyverse)
library(here)
library(loo)

# setwd("/u/home/t/tagibson/Projects/CovSelectMeta")
# i_am("AB-NMA-dong.R")

rstan_options("auto_write" = TRUE)
options(mc.cores = 4)

# compile models

# nma.reglambda.sepvar <- stan_model(here("Models", "AB-NMA-sep-var-reglambda.stan"))
nma.reglambda.sepvar.freegamma <- stan_model(here("Models", "AB-NMA-sepvar-hier-freegamma.stan"))

#######
########

set.seed(714) # like babe ruth
nma.dat.1 <- read_csv(here("Data", "dat.dong2013.csv")) # %>%
  # group_by(treatment) %>%
  # filter(n() > 2) %>%
  # ungroup() %>% group_by(id) %>%
  # filter(n() > 1)

len <- dim(nma.dat.1)[1]
ntrt <- length(unique(nma.dat.1$treatment))
nstudy <- max(nma.dat.1$id)

nma.dat.1$treatment <- factor(nma.dat.1$treatment,
                              levels = names(table(nma.dat.1$treatment))[order(table(nma.dat.1$treatment))])

t <- as.numeric(nma.dat.1$treatment)
s <- nma.dat.1$id
r <- nma.dat.1$death
totaln <- nma.dat.1$randomized
zeros <- rep(0, ntrt)
Lambda <- diag(ntrt)
nu <- ntrt + 1

######## raw summary of data
# nma.dat.1 %>%
#   group_by(treatment) %>%
#   mutate(prob = (death + 0.5) / (randomized + 0.5),
#          logit.prob = log(prob / (1 - prob)),
#          mu.t = mean(logit.prob)) %>%
#   summarize(mu.t = mean(mu.t),
#             studies = n())


tau0_sepvar <- unlist(nma.dat.1 %>%
                        group_by(treatment) %>%
                        mutate(prob = (death + 0.5) / (randomized + 0.5),
                               logit.prob = log(prob / (1 - prob)),
                               mu.t = mean(logit.prob)) %>%
                        ungroup() %>%
                        group_by(id) %>%
                        mutate(beta.i = mean(logit.prob - mu.t)) %>%
                        ungroup() %>%
                        mutate(nu.it = logit.prob - (mu.t + beta.i)) %>%
                        group_by(treatment) %>%
                        summarize(sd.nu = sd(nu.it),
                                  sd.sd.nu = sd.nu * gamma((n() - 1) / 2) / gamma(n() / 2) * sqrt((n() - 1) / 2 - (gamma(n() / 2) / gamma((n() - 1) / 2))^2)) %>%
                        summarize(tau0 = (n() - 2) / 2 * mean(sd.sd.nu)))

reglambda.sepvar.freegamma.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                                       s = s, r = r, totaln = totaln,
                                       higher_better = 0, zeros = zeros,
                                       scale_gamma = 1,
                                       scale_global_lambda = tau0_sepvar, nu_global_lambda = 3, nu_local_lambda = 1,
                                       slab_scale_lambda = 2, slab_df_lambda = 4)




# fit.1.reglambda.sepvar <- sampling(nma.reglambda.sepvar,
#                                    pars = c("SD", "sigma_beta", "CORR",
#                                             "MU", "AR", "LOR", "log_likelihood"),
#                                    data = reglambda.dat, 
#                                    iter = 6000, warmup = 2000, thin = 2, cores = 4,
#                                    #iter = 400, warmup = 200, cores = 4,
#                                    control = list(adapt_delta = 0.99, max_treedepth = 15),
#                                    seed = 923)

RHS.NS.SV <- sampling(nma.reglambda.sepvar.freegamma,
                                   pars = c("SD", "sigma_beta", "CORR",
                                            "MU", "AR", "LOR", "log_likelihood"),
                                   data = reglambda.sepvar.freegamma.dat,
                                   iter = 6000, warmup = 2000, thin = 2, cores = 4,
                                   #iter = 100, warmup = 50, cores = 4,
                                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                                   seed = 919)



round(summary(RHS.NS.SV, pars =      c("SD", "CORR", "lp__"))$summary, 4)


log_lik <-  extract_log_lik(RHS.NS.SV, parameter_name = "log_likelihood")

r_eff <- relative_eff(exp(log_lik), chain_id = rep(1:4, each = 2000))
# r_eff <- lapply(log_lik, function(x){
#   relative_eff(exp(x), chain_id = rep(1:4, each = 3000)) # 4 chains, each with 3000 draws
# })

loo_1 <- loo(log_lik, r_eff = r_eff)
# loo_list <- lapply(1:length(log_lik), function(j){
#   loo::loo(log_lik[[j]], r_eff = r_eff[[j]])
# })


#names(models) <- names(loo_list) <- c("sepvar.freegamma")

write_rds(list(RHS.NS.SV, loo_1), file = here("Results", "dong-models6-allarms.rds"))
