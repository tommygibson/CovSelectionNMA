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
nma.std <- stan_model(here("Models", "AB-NMA.stan"))
nma.lkj <- stan_model(here("Models", "AB-NMA-LKJ.stan"))
nma.reglambda <- stan_model(here("Models", "AB-NMA-regularized-lambda.stan"))
nma.RHS.NS <- stan_model(here("Models", "AB-NMA-RHS-NS.stan"))
nma.reglambda.sepvar <- stan_model(here("Models", "AB-NMA-sep-var-reglambda.stan"))
nma.reglambda.sepvar.hier <- stan_model(here("Models", "AB-NMA-sep-var-hierarchical.stan"))
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

tau0 <- unlist(nma.dat.1 %>%
  group_by(treatment) %>%
  mutate(logit.prob = log(((death + 0.5) / (randomized + 0.5)) / (((randomized + 0.5) - (death + 0.5)) / (randomized + 0.5))),
         weight = randomized) %>%
  summarize(logit.mean = mean(logit.prob),
            logit.mean.weighted = weighted.mean(logit.prob, w = weight),
            logit.sd = sd(logit.prob),
            logit.sd.sd = logit.sd * gamma((n() - 1) / 2) / gamma(n() / 2) * sqrt((n() - 1) / 2 - (gamma(n() / 2) / gamma((n() - 1) / 2))^2),
            studies = n()) %>%
  summarize(tau0 = (n() - 2) / 2 * mean(logit.sd.sd, na.rm = TRUE)))


std.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                   s = s, r = r, totaln = totaln,
                   higher_better = 0, zeros = zeros,
                   Lambda = Lambda, 
                   nu = nu)
lkj.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                     s = s, r = r, totaln = totaln,
                     higher_better = 0, zeros = zeros,
                     Lambda = Lambda, 
                     nu = nu, eta = 1)

RHS.CS.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                      s = s, r = r, totaln = totaln,
                      higher_better = 0, zeros = zeros,
                      scale_gamma = sqrt(2),
                      scale_global_lambda = tau0, nu_global_lambda = 3, nu_local_lambda = 1,
                      slab_scale_lambda = 2, slab_df_lambda = 4)
RHS.NS.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                      s = s, r = r, totaln = totaln,
                      higher_better = 0, zeros = zeros,
                      scale_gamma = 1,
                      scale_global_lambda = tau0, nu_global_lambda = 3, nu_local_lambda = 1,
                      slab_scale_lambda = 2, slab_df_lambda = 4)
RHS.CS.SV.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                             s = s, r = r, totaln = totaln,
                             higher_better = 0, zeros = zeros,
                             scale_gamma = 4,
                             scale_global_lambda = tau0_sepvar, nu_global_lambda = 3, nu_local_lambda = 1,
                             slab_scale_lambda = 2, slab_df_lambda = 4)
RHS.NS.SV.dat <- list(len = len, ntrt = ntrt, nstudy = nstudy, t = t,
                                      s = s, r = r, totaln = totaln,
                                      higher_better = 0, zeros = zeros,
                                      scale_gamma = 1,
                                      scale_global_lambda = tau0_sepvar, nu_global_lambda = 3, nu_local_lambda = 1,
                                      slab_scale_lambda = 2, slab_df_lambda = 4)


# 
# fit.1.iw <- sampling(nma.std, pars = c("SD", "CORR", "MU", "AR",
#                                          "LOR", "log_likelihood"), 
#                        data = std.dat, iter = 5000, warmup = 2000, thin = 1, cores = 4,
#                        #data = std.dat, iter = 400, warmup = 200, cores = 4,
#                        control = list(adapt_delta = 0.99),
#                        seed = 919)
# 
# fit.1.lkj <- sampling(nma.lkj, pars = c("SD", "CORR", "MU", "AR",
#                                         "LOR", "log_likelihood"), 
#                       data = lkj.dat, iter = 5000, warmup = 2000, thin = 1, cores = 4,
#                       #data = lkj.dat, iter = 400, warmup = 200, cores = 4,
#                       control = list(adapt_delta = 0.99, max_treedepth = 15),
#                       seed = 920)
# 
# fit.1.reglambda <- sampling(nma.reglambda, 
#                             pars = c("SD", "CORR", "MU", "AR",
#                                      "LOR", "log_likelihood"), 
#                             data = reglambda.dat, iter = 5000, warmup = 2000, thin = 1, cores = 4,
#                             #data = reglambda.dat, iter = 400, warmup = 200, cores = 4,
#                             control = list(adapt_delta = 0.99, max_treedepth = 15),
#                             seed = 921)
# fit.1.RHS.NS <- sampling(nma.RHS.NS, 
#                          pars = c("SD", "CORR", "MU", "AR",
#                                   "LOR", "log_likelihood"), 
#                          data = RHS.NS.dat, iter = 5000, warmup = 2000, thin = 1, cores = 4,
#                          #data = RHS.NS.dat, iter = 400, warmup = 200, cores = 4,
#                          control = list(adapt_delta = 0.99, max_treedepth = 15),
#                          seed = 922)
# fit.1.reglambda.sepvar <- sampling(nma.reglambda.sepvar,
#                                    pars = c("SD", "sigma_beta", "CORR",
#                                             "MU", "AR", "LOR", "log_likelihood"),
#                                    data = reglambda.dat, 
#                                    iter = 6000, warmup = 2000, thin = 1, cores = 4,
#                                    #iter = 400, warmup = 200, cores = 4,
#                                    control = list(adapt_delta = 0.99, max_treedepth = 15),
#                                    seed = 923)

# fit.1.sepvar.freegamma <- sampling(nma.reglambda.sepvar.hier.freegamma,
#                               pars = c("SD", "sigma_beta", "CORR",
#                                        "MU", "AR", "LOR", "log_likelihood"),
#                               data = reglambda.sepvar.freegamma.dat,
#                               iter = 6000, warmup = 2000, thin = 2, cores = 4,
#                               #iter = 400, warmup = 200, cores = 4,
#                               control = list(adapt_delta = 0.99, max_treedepth = 15),
#                               seed = 924)

iw.ordered <- sampling(nma.std, pars = c("SD", "CORR", "MU", "AR",
                                         "LOR", "log_likelihood"), 
                       data = std.dat, iter = 6000, warmup = 2000, thin = 2, cores = 4,
                       #data = std.dat, iter = 400, warmup = 200, cores = 4,
                       control = list(adapt_delta = 0.99),
                       seed = 919)

lkj.ordered <- sampling(nma.lkj, pars = c("SD", "CORR", "MU", "AR",
                                          "LOR", "log_likelihood"), 
                        data = lkj.dat, iter = 6000, warmup = 2000, thin = 2, cores = 4,
                        #data = lkj.dat, iter = 400, warmup = 200, cores = 4,
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        seed = 919)

RHS.CS.ordered <- sampling(nma.reglambda, 
                        pars = c("SD", "CORR", "MU", "AR",
                                 "LOR", "log_likelihood"), 
                        data = RHS.CS.dat, iter = 6000, warmup = 2000, thin = 2, cores = 4,
                        #data = reglambda.dat, iter = 400, warmup = 200, cores = 4,
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        seed = 919)
RHS.NS.ordered <- sampling(nma.RHS.NS, 
                           pars = c("SD", "CORR", "MU", "AR",
                                    "LOR", "log_likelihood"), 
                           data = RHS.NS.dat, iter = 6000, warmup = 2000, thin = 2, cores = 4,
                           #data = RHS.NS.dat, iter = 400, warmup = 200, cores = 4,
                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                           seed = 919)
RHS.CS.SV.ordered <- sampling(nma.reglambda.sepvar,
                           pars = c("SD", "sigma_beta", "CORR",
                                    "MU", "AR", "LOR", "log_likelihood"),
                           data = RHS.CS.SV.dat, 
                           iter = 6000, warmup = 2000, thin = 2, cores = 4,
                           #iter = 400, warmup = 200, cores = 4,
                           control = list(adapt_delta = 0.99, max_treedepth = 15),
                           seed = 919)


round(summary(iw.ordered, pars =      c("SD", "CORR", "lp__"))$summary, 4)
round(summary(lkj.ordered, pars =       c("SD", "CORR", "lp__"))$summary, 4)
round(summary(RHS.CS.ordered, pars =     c("SD", "CORR", "lp__"))$summary, 4)
round(summary(RHS.NS.ordered, pars =  c("SD", "CORR", "lp__"))$summary, 4)
round(summary(RHS.CS.SV.ordered, pars =  c("SD", "sigma_beta", "CORR", "lp__"))$summary, 4)
round(summary(iw.ordered, pars =      c("MU", "AR", "lp__"))$summary, 4)
round(summary(lkj.ordered, pars =       c("MU", "AR", "lp__"))$summary, 4)
round(summary(RHS.CS.ordered, pars =     c("MU", "AR", "lp__"))$summary, 4)
round(summary(RHS.NS.ordered, pars =  c("MU", "AR", "lp__"))$summary, 4)
round(summary(RHS.CS.SV.ordered, pars =  c("MU", "AR", "lp__"))$summary, 4)

log_lik <- list()
log_lik[[1]] <- extract_log_lik(iw.ordered, parameter_name = "log_likelihood")
log_lik[[2]] <- extract_log_lik(lkj.ordered, parameter_name = "log_likelihood")
log_lik[[3]] <- extract_log_lik(RHS.CS.ordered, parameter_name = "log_likelihood")
log_lik[[4]] <- extract_log_lik(RHS.NS.ordered, parameter_name = "log_likelihood")
log_lik[[5]] <- extract_log_lik(RHS.CS.SV.ordered, parameter_name = "log_likelihood")
#log_lik[[5]] <- extract_log_lik(fit.1.sepvar.freegamma, parameter_name = "log_likelihood")

r_eff <- lapply(log_lik, function(x){
  relative_eff(exp(x), chain_id = rep(1:4, each = 4000)) # 4 chains, each with 3000 draws
})

loo_list <- lapply(1:length(log_lik), function(j){
  loo::loo(log_lik[[j]], r_eff = r_eff[[j]])
})

models <- list(fit.1.stan, fit.1.lkj, fit.1.reglambda, fit.1.RHS.NS, fit.1.reglambda.sepvar) # fit.1.sepvar.freegamma)

names(models) <- names(loo_list) <- c("IW", "LKJ", "RHS-CS", "RHS-NS", "RHS-CS-SV")

write_rds(list(models, loo_list), file = here("Results", "dong-models1-5-allarms.rds"))
