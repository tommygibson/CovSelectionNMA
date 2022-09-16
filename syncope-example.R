##### 3RE model with covariance selection for syncope data

library(rstan)
library(tidyverse)
library(here)
library(loo)
library(extrafont)
library(xtable)

source(here("CTS-functions.R"))

options(mc.cores = parallel::detectCores())
# define the models

re3.std.model <- stan_model(here("Models", "3RE.stan"))
re3.iw.model <- stan_model(here("Models", "3RE-IW.stan"))
re3.lkj.model <- stan_model(here("Models", "3RE-LKJ.stan"))
re3.reglambda.model <- stan_model(here("Models", "3RE-regularized-lambda-cholesky.stan"))
re3.freegamma.model <- stan_model(here("Models", "3RE-reglambda-freegamma.stan"))

# load syncope data, filter out vars with only 2-4 observations
# also remove studies without counts

syncope <- read_csv(here("Data", "syncope.cleaned.csv")) %>%
  group_by(Variable) %>%
  filter(!is.na(n_i0)) %>%
  filter(n() > 2)

syncope_re_hats <- syncope %>%
  transmute(Variable = Variable,
            nu.hat = log((n_i1 / N_i) / (1 - n_i1 / N_i)),
            beta.hat = 1 / 2 * (log((y_i1 / n_i1) / (1 - y_i1 / n_i1)) + log((y_i0 / n_i0) / (1 - y_i0 / n_i0))),
            delta.hat = lnORhat) %>%
  group_by(Variable) %>%
  summarize(sd.beta.hat = sd(beta.hat),
            sd.sd.beta.hat = sd.beta.hat * sqrt(1 - 2 / (n() - 1) * (gamma(n() / 2) / gamma((n() - 1) / 2))^2),
            sd.delta.hat = sd(delta.hat),
            sd.sd.delta.hat = sd.delta.hat * sqrt(1 - 2 / (n() - 1) * (gamma(n() / 2) / gamma((n() - 1) / 2))^2),
            sd.nu.hat = sd(nu.hat),
            sd.sd.nu.hat = sd.nu.hat * sqrt(1 - 2 / (n() - 1) * (gamma(n() / 2) / gamma((n() - 1) / 2))^2),
            tau0 = 2 * mean(c(sd.sd.beta.hat, sd.sd.nu.hat, sd.sd.delta.hat))) %>%
  filter(sd.nu.hat < 0.1 | sd.beta.hat < 0.1 | sd.delta.hat < 0.1) 

# print(xtable(syncope_re_hats, type = "latex"),
#       file = here("Results", "syncope-small-sds.tex"), include.rownames = FALSE)

# subset syncope data to only those with small RESDs
syncope_small_sds <- syncope %>%
  filter(Variable %in% syncope_re_hats$Variable)
  
  small_sd <- syncope_small_sds %>%
    ungroup %>%
    filter(Variable == "Chest Pain")
  
  y <- as.matrix(small_sd %>% 
    dplyr::select(y_i1, y_i0) )
  n <- as.matrix(small_sd %>% 
    dplyr::select(n_i1, n_i0))
  bigN <- rowSums(n)
  S <- nrow(y)
  
  nu.hat <- log((n[,1] / bigN) / (1 - n[,1] / bigN))
  beta.hat <- apply(log((y/n) / (1 - y/n)), 1, mean)
  delta.hat <- small_sd$lnORhat
  
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
  max.sd.sd <- max(c(sd.sd.nu.hat, sd.sd.beta.hat, sd.sd.delta.hat))
  scale_global_lambda <- 2 / (3 - 2) * mean.sd.sd
  scale_global_lambda_max <- 2 / (3 - 2) * max.sd.sd
  
  re3.dat.stan <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                       a = 0, b = 2, c = 0, d = 2, f = 0, g = 2)
  re3.dat.lkj <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                      a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                      nu = 1)
  re3.dat.reglambda <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                           a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                           scale_global_lambda = scale_global_lambda, scale_gamma = 2,
                           nu_global_lambda = 3,  nu_local_lambda = 1,
                           slab_scale_lambda = 1, slab_df_lambda = 4)
  re3.dat.freegamma <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                            a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                            scale_global_lambda = scale_global_lambda, scale_gamma = 1 / sqrt(2),
                            nu_global_lambda = 3,  nu_local_lambda = 1,
                            slab_scale_lambda = 1, slab_df_lambda = 4)
  # scale_global_lambda = 1, scale_global_gamma = 1,
  # nu_global_lambda = 1, nu_global_gamma = 1,
  # nu_local_lambda = 1, nu_local_gamma = 1,
  # slab_scale_lambda = 2, slab_df_lambda = 4,
  # slab_scale_gamma = 1, slab_df_gamma = 8
  re3.stan <- sampling(re3.std.model, data = re3.dat.stan,
                       chains = 4, iter = 4000,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       seed = 714)
  re3.iw.stan <- sampling(re3.iw.model, data = re3.dat.stan,
                          chains = 4, iter = 4000,
                          pars = c("theta0", "SD", "CORR", "Sigma"),
                          seed = 715)#,
                          #control = list(adapt_delta = 0.95))
  
  re3.lkj.stan <- sampling(re3.lkj.model, data = re3.dat.lkj,
                           chains = 4, iter = 4000, warmup = 2000,
                           pars = c("theta0", "SD", "CORR", "Sigma"),
                           control = list(adapt_delta = .99, max_treedepth = 15),
                           cores = 4,
                           seed = 716)
  
  re3.reglambda <- sampling(re3.reglambda.model, data = re3.dat.reglambda,
                            chains = 4, iter = 4000, warmup = 2000,
                            pars = c("theta0", "SD", "CORR", "Sigma"),
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            cores = 4,
                            seed = 717)
  re3.freegamma <- sampling(re3.freegamma.model, data = re3.dat.freegamma,
                            chains = 4, iter = 4000, warmup = 2000,
                            pars = c("theta0", "SD", "CORR", "Sigma"),
                            control = list(adapt_delta = 0.99, max_treedepth = 15),
                            cores = 4,
                            seed = 718)
  
  round(summary(re3.stan, pars = c("beta0", "delta0", "nu0", "sigma_beta", "sigma_delta", "sigma_nu", "lp__"))$summary, 4)
  round(summary(re3.iw.stan, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
  round(summary(re3.lkj.stan, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
  round(summary(re3.reglambda, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
  round(summary(re3.freegamma, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
  

sd.sims <- cbind.data.frame(c(#as.vector(do.call(cbind, rstan::extract(re3.stan, pars = c("sigma_beta", "sigma_delta", "sigma_nu")))),
                              as.vector(rstan::extract(re3.iw.stan, pars = "SD")$SD),
                              as.vector(rstan::extract(re3.lkj.stan, pars = "SD")$SD),
                              as.vector(rstan::extract(re3.reglambda, pars = "SD")$SD)),
                            rep(rep(c("sigma[beta]", "sigma[delta]", "sigma[nu]"), each = 8000), 3),
                            rep(c("IW", "LKJ", "RHS-CS"), each = 24000))
names(sd.sims) <- c("samples", "parameter", "model")

# extract only unique correlation values rather than whole 3x3 matrix
# only need rho_beta_delta, rho_beta_nu, rho_delta_nu, which are in the 2nd, 3rd, and 6th spots

corr.sims <- cbind.data.frame(c(as.vector(matrix(rstan::extract(re3.iw.stan, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)]),
                                as.vector(matrix(rstan::extract(re3.lkj.stan, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)]),
                                as.vector(matrix(rstan::extract(re3.reglambda, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)])),
                              rep(rep(c("rho[beta][delta]", "rho[beta][nu]", "rho[delta][nu]"), each = 8000), 3),
                              rep(c("IW", "LKJ", "RHS-CS"), each = 24000))
names(corr.sims) <- c("samples", "parameter", "model")


sd.means <- sd.sims %>%
  ungroup() %>%
  group_by(parameter, model) %>%
  summarize(mean = round(mean(samples), 3),
            median = round(quantile(samples, .5), 3),
            x.mean = 1,
            y.mean = max(ifelse(str_detect(parameter, "beta"), 775, 
                                ifelse(str_detect(parameter, "delta"), 850, 2100))),
            x.median = 1,
            y.median = max(ifelse(str_detect(parameter, "beta"), 700, 
                                  ifelse(str_detect(parameter, "delta"), 775, 1900))))

corr.means <- corr.sims %>%
  ungroup() %>%
  group_by(parameter, model) %>%
  summarize(mean = round(mean(samples), 3),
            sd = round(sd(samples), 3),
            x.mean = -1,
            y.mean = max(ifelse(str_detect(parameter, "delta"), 1750, 2100)),
            x.sd = -1,
            y.sd = max(ifelse(str_detect(parameter, "delta"), 1600, 1900)))

sd.chest.plot <- sd.sims %>%
  ggplot(aes(x = samples)) +
  geom_histogram(bins = 30) +
  facet_grid(parameter ~ model,
             scales = "free",
             labeller = label_parsed) +
  xlim(c(0, 2)) +
  theme_bw() +
  theme(text = element_text(family =  "LM Roman 10"))  +
  geom_text(data = sd.means, aes(x = x.mean, y = y.mean, label = paste0("Mean = ", mean)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  geom_text(data = sd.means, aes(x = x.median, y = y.median, label = paste0("Med. = ", median)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  xlab("SD") +
  ylab("Samples") 
  

corr.chest.plot <- corr.sims %>%
  ggplot(aes(x = samples)) + 
  geom_histogram(bins = 40) + 
  facet_grid(parameter ~ model,
             scales = "free",
             labeller = label_parsed) +
  xlim(c(-1, 1)) +
  theme_bw() +
  theme(text = element_text(family = "LM Roman 10")) +
  geom_text(data = corr.means, aes(x = x.mean, y = y.mean, label = paste0("Mean = ", mean)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  geom_text(data = corr.means, aes(x = x.sd, y = y.sd, label = paste0("SD = ", sd)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  xlab("Correlation") +
  ylab("Samples")

# ggsave(here("Figures", "sd_chest_plot.pdf"), plot = sd.chest.plot,
#        height = 5, width = 5, 
#        units = "in", device = "pdf", dpi = 300)
# ggsave(here("Figures", "corr_chest_plot.pdf"), plot = corr.chest.plot,
#        height = 5, width = 5,
#        units = "in", device = "pdf", dpi = 300)
    

reglambda.sims <- make_mean_var_sims(rstan::extract(re3.reglambda, pars = c("theta0", "Sigma")))
lkj.sims <- make_mean_var_sims(rstan::extract(re3.lkj.stan, pars = c("theta0", "Sigma")))
iw.sims <- make_mean_var_sims(rstan::extract(re3.iw.stan, pars = c("theta0", "Sigma")))
cts0_from_one_posterior(reglambda.sims[1,])
reglambda.cts0.summary <- as.data.frame(cts_summary_from_gamma(reglambda.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS =  V1,
         Mean = V2,
         SD =   V3)
lkj.cts0.summary <- as.data.frame(cts_summary_from_gamma(lkj.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS = V1,
         Mean = V2,
         SD = V3)
iw.cts0.summary <- as.data.frame(cts_summary_from_gamma(iw.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS = V1,
         Mean = V2,
         SD = V3)


reglambda.cts0 <- reglambda.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
lkj.cts0 <- lkj.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
iw.cts0 <- iw.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
all.cts0 <- cbind.data.frame(iw.cts0, lkj.cts0$Posterior, reglambda.cts0$Posterior)
names(all.cts0) <- c("CTS", "IW", "LKJ", "RHS-CS")

# print(xtable(all.cts0, type = "latex"),
#       file = here("Results", "chest-pain-results.tex"), include.rownames = FALSE)


########################################################################
######## MALE GENDER
########################################################################
small_sd <- syncope_small_sds %>%
  ungroup %>%
  filter(Variable == "Male Gender")

y <- as.matrix(small_sd %>% 
                 dplyr::select(y_i1, y_i0) )
n <- as.matrix(small_sd %>% 
                 dplyr::select(n_i1, n_i0))
bigN <- rowSums(n)
S <- nrow(y)

nu.hat <- log((n[,1] / bigN) / (1 - n[,1] / bigN))
beta.hat <- apply(log((y/n) / (1 - y/n)), 1, mean)
delta.hat <- small_sd$lnORhat

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
max.sd.sd <- max(c(sd.sd.nu.hat, sd.sd.beta.hat, sd.sd.delta.hat))
scale_global_lambda <- 2 / (3 - 2) * mean.sd.sd
scale_global_lambda_max <- 2 / (3 - 2) * max.sd.sd

re3.dat.stan <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                     a = 0, b = 2, c = 0, d = 2, f = 0, g = 2)
re3.dat.lkj <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                    a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                    nu = 1)
re3.dat.reglambda <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                          a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                          scale_global_lambda = scale_global_lambda, scale_gamma = 4,
                          nu_global_lambda = 3,  nu_local_lambda = 1,
                          slab_scale_lambda = 1, slab_df_lambda = 4)
re3.dat.freegamma <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                          a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                          scale_global_lambda = scale_global_lambda, scale_gamma = 1 / sqrt(2),
                          nu_global_lambda = 3,  nu_local_lambda = 1,
                          slab_scale_lambda = 1, slab_df_lambda = 4)

re3.stan <- sampling(re3.std.model, data = re3.dat.stan,
                     chains = 4, iter = 4000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     seed = 714)

re3.iw.stan <- sampling(re3.iw.model, data = re3.dat.stan,
                        chains = 4, iter = 4000,
                        pars = c("theta0", "SD", "CORR", "Sigma"),
                        seed = 715)

re3.lkj.stan <- sampling(re3.lkj.model, data = re3.dat.lkj,
                         chains = 4, iter = 6000, warmup = 4000,
                         pars = c("theta0", "SD", "CORR", "Sigma"),
                         control = list(adapt_delta = .99, max_treedepth = 15),
                         cores = 4,
                         seed = 716)

re3.reglambda <- sampling(re3.reglambda.model, data = re3.dat.reglambda,
                          chains = 4, iter = 5000, warmup = 3000,
                          pars = c("theta0", "SD", "CORR", "Sigma"),
                          control = list(adapt_delta = 0.99, max_treedepth = 15),
                          cores = 4,
                          seed = 717)

re3.freegamma <- sampling(re3.freegamma.model, data = re3.dat.freegamma,
                          chains = 4, iter = 4000, warmup = 2000,
                          pars = c("theta0", "SD", "CORR", "Sigma"),
                          control = list(adapt_delta = 0.99, max_treedepth = 15),
                          cores = 4,
                          seed = 718)

round(summary(re3.stan, pars = c("beta0", "delta0", "nu0", "sigma_beta", "sigma_delta", "sigma_nu", "lp__"))$summary, 4)
round(summary(re3.iw.stan, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
round(summary(re3.lkj.stan, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
round(summary(re3.reglambda, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
round(summary(re3.freegamma, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)


sd.sims <- cbind.data.frame(c(as.vector(rstan::extract(re3.iw.stan, pars = "SD")$SD),
                              as.vector(rstan::extract(re3.lkj.stan, pars = "SD")$SD),
                              as.vector(rstan::extract(re3.reglambda, pars = "SD")$SD)),
                            rep(rep(c("sigma[beta]", "sigma[delta]", "sigma[nu]"), each = 8000), 3),
                            rep(c("IW", "LKJ", "RHS-CS"), each = 24000))
names(sd.sims) <- c("samples", "parameter", "model")

# extract only unique correlation values rather than whole 3x3 matrix
# only need rho_beta_delta, rho_beta_nu, rho_delta_nu, which are in the 2nd, 3rd, and 6th spots

corr.sims <- cbind.data.frame(c(as.vector(matrix(rstan::extract(re3.iw.stan, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)]),
                                as.vector(matrix(rstan::extract(re3.lkj.stan, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)]),
                                as.vector(matrix(rstan::extract(re3.reglambda, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)])),
                              rep(rep(c("rho[beta][delta]", "rho[beta][nu]", "rho[delta][nu]"), each = 8000), 3),
                              rep(c("IW", "LKJ", "RHS-CS"), each = 24000))
names(corr.sims) <- c("samples", "parameter", "model")


sd.means <- sd.sims %>%
  ungroup() %>%
  group_by(parameter, model) %>%
  summarize(mean = round(mean(samples), 3),
            median = round(quantile(samples, .5), 3),
            x.mean = 1,
            y.mean = max(ifelse(str_detect(parameter, "beta"), 775, 
                                ifelse(str_detect(parameter, "delta"), 2350, 4700))),
            x.median = 1,
            y.median = max(ifelse(str_detect(parameter, "beta"), 725, 
                                  ifelse(str_detect(parameter, "delta"), 2150, 4300))))

corr.means <- corr.sims %>%
  ungroup() %>%
  group_by(parameter, model) %>%
  summarize(mean = round(mean(samples), 3),
            sd = round(sd(samples), 3),
            x.mean = -1,
            y.mean = max(ifelse(str_detect(parameter, "beta"), 1425, 2100)),
            x.sd = -1,
            y.sd = max(ifelse(str_detect(parameter, "beta"), 1325, 1900)))

sd.male.plot <- sd.sims %>%
  ggplot(aes(x = samples)) +
  geom_histogram(bins = 40) +
  facet_grid(parameter ~ model,
             scales = "free",
             labeller = label_parsed) +
  xlim(c(0, 2)) +
  theme_bw() +
  theme(text = element_text(family =  "LM Roman 10")) +
  geom_text(data = sd.means, aes(x = x.mean, y = y.mean, label = paste0("Mean = ", mean)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  geom_text(data = sd.means, aes(x = x.median, y = y.median, label = paste0("Med. = ", median)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  xlab("SD") +
  ylab("Samples")


corr.male.plot <- corr.sims %>%
  ggplot(aes(x = samples)) + 
  geom_histogram(bins = 40) + 
  facet_grid(parameter ~ model,
             scales = "free",
             labeller = label_parsed) +
  xlim(c(-1, 1)) +
  theme_bw() +
  theme(text = element_text(family = "LM Roman 10")) +
  geom_text(data = corr.means, aes(x = x.mean, y = y.mean, label = paste0("Mean = ", mean)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  geom_text(data = corr.means, aes(x = x.sd, y = y.sd, label = paste0("SD = ", sd)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  xlab("Correlation") +
  ylab("Samples")

# ggsave(here("Figures", "sd_male_plot.pdf"), plot = sd.male.plot,
#        height = 5, width = 5, units = "in", device = "pdf", dpi = 300)
# ggsave(here("Figures", "corr_male_plot.pdf"), plot = corr.male.plot,
#        height = 5, width = 5, units = "in", device = "pdf", dpi = 300)

reglambda.sims <- make_mean_var_sims(rstan::extract(re3.reglambda, pars = c("theta0", "Sigma")))
freegamma.sims <- make_mean_var_sims(rstan::extract(re3.freegamma, pars = c("theta0", "Sigma")))
lkj.sims <- make_mean_var_sims(rstan::extract(re3.lkj.stan, pars = c("theta0", "Sigma")))
iw.sims <- make_mean_var_sims(rstan::extract(re3.iw.stan, pars = c("theta0", "Sigma")))
cts0_from_one_posterior(reglambda.sims[1,])
reglambda.cts0.summary <- as.data.frame(cts_summary_from_gamma(reglambda.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS =  V1,
         Mean = V2,
         SD =   V3)
freegamma.cts0.summary <- as.data.frame(cts_summary_from_gamma(freegamma.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS =  V1,
         Mean = V2,
         SD =   V3)
lkj.cts0.summary <- as.data.frame(cts_summary_from_gamma(lkj.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS = V1,
         Mean = V2,
         SD = V3)
iw.cts0.summary <- as.data.frame(cts_summary_from_gamma(iw.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS = V1,
         Mean = V2,
         SD = V3)

reglambda.cts0 <- reglambda.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
freegamma.cts0 <- freegamma.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
lkj.cts0 <- lkj.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
iw.cts0 <- iw.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
all.cts0 <- cbind.data.frame(iw.cts0, lkj.cts0$Posterior, reglambda.cts0$Posterior)
names(all.cts0) <- c("CTS", "IW", "LKJ", "RHS-CS")

# print(xtable(all.cts0, type = "latex"),
#       file = here("Results", "male-gender-results.tex"), include.rownames = FALSE)


########################################################################
######## WHITE RACE
########################################################################

small_sd <- syncope_small_sds %>%
  ungroup %>%
  filter(Variable == "White Race")

y <- as.matrix(small_sd %>% 
                 dplyr::select(y_i1, y_i0) )
n <- as.matrix(small_sd %>% 
                 dplyr::select(n_i1, n_i0))
bigN <- rowSums(n)
S <- nrow(y)

nu.hat <- log((n[,1] / bigN) / (1 - n[,1] / bigN))
beta.hat <- apply(log((y/n) / (1 - y/n)), 1, mean)
delta.hat <- small_sd$lnORhat

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
max.sd.sd <- max(c(sd.sd.nu.hat, sd.sd.beta.hat, sd.sd.delta.hat))
scale_global_lambda <- 2 / (3 - 2) * mean.sd.sd
scale_global_lambda_max <- 2 / (3 - 2) * max.sd.sd

re3.dat.stan <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                     a = 0, b = 2, c = 0, d = 2, f = 0, g = 2)
re3.dat.lkj <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                    a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                    nu = 1)
re3.dat.reglambda <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                          a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                          scale_global_lambda = scale_global_lambda, scale_gamma = 2,
                          nu_global_lambda = 3,  nu_local_lambda = 1,
                          slab_scale_lambda = 1, slab_df_lambda = 4)
re3.dat.freegamma <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                          a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                          scale_global_lambda = scale_global_lambda, scale_gamma = 1 / sqrt(2),
                          nu_global_lambda = 3,  nu_local_lambda = 1,
                          slab_scale_lambda = 1, slab_df_lambda = 4)

re3.stan <- sampling(re3.std.model, data = re3.dat.stan,
                     chains = 4, iter = 4000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     seed = 714)
re3.iw.stan <- sampling(re3.iw.model, data = re3.dat.stan,
                        chains = 4, iter = 4000,
                        pars = c("theta0", "SD", "CORR", "Sigma"),
                        seed = 715)

re3.lkj.stan <- sampling(re3.lkj.model, data = re3.dat.lkj,
                         chains = 4, iter = 6000, warmup = 4000,
                         pars = c("theta0", "SD", "CORR", "Sigma"),
                         control = list(adapt_delta = .995, max_treedepth = 15),
                         cores = 4,
                         seed = 716)

re3.reglambda <- sampling(re3.reglambda.model, data = re3.dat.reglambda,
                          chains = 4, iter = 6000, warmup = 4000,
                          pars = c("theta0", "SD", "CORR", "Sigma"),
                          control = list(adapt_delta = 0.99, max_treedepth = 15),
                          cores = 4,
                          seed = 717)
re3.freegamma <- sampling(re3.freegamma.model, data = re3.dat.freegamma,
                          chains = 4, iter = 6000, warmup = 4000,
                          pars = c("theta0", "SD", "CORR", "Sigma"),
                          control = list(adapt_delta = 0.995, max_treedepth = 15),
                          cores = 4,
                          seed = 718)

round(summary(re3.stan, pars = c("beta0", "delta0", "nu0", "sigma_beta", "sigma_delta", "sigma_nu", "lp__"))$summary, 4)
round(summary(re3.iw.stan, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
round(summary(re3.lkj.stan, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
round(summary(re3.reglambda, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)
round(summary(re3.freegamma, pars = c("theta0", "SD", "CORR", "lp__"))$summary, 4)


sd.sims <- cbind.data.frame(c(#as.vector(do.call(cbind, rstan::extract(re3.stan, pars = c("sigma_beta", "sigma_delta", "sigma_nu")))),
  as.vector(rstan::extract(re3.iw.stan, pars = "SD")$SD),
  as.vector(rstan::extract(re3.lkj.stan, pars = "SD")$SD),
  as.vector(rstan::extract(re3.reglambda, pars = "SD")$SD)),
  rep(rep(c("sigma[beta]", "sigma[delta]", "sigma[nu]"), each = 8000), 3),
  rep(c("IW", "LKJ", "RHS-CS"), each = 24000))
names(sd.sims) <- c("samples", "parameter", "model")

# extract only unique correlation values rather than whole 3x3 matrix
# only need rho_beta_delta, rho_beta_nu, rho_delta_nu, which are in the 2nd, 3rd, and 6th spots

corr.sims <- cbind.data.frame(c(as.vector(matrix(rstan::extract(re3.iw.stan, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)]),
                                as.vector(matrix(rstan::extract(re3.lkj.stan, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)]),
                                as.vector(matrix(rstan::extract(re3.reglambda, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)])),
                              rep(rep(c("rho[beta][delta]", "rho[beta][nu]", "rho[delta][nu]"), each = 8000), 3),
                              rep(c("IW", "LKJ", "RHS-CS"), each = 24000))
names(corr.sims) <- c("samples", "parameter", "model")


sd.means <- sd.sims %>%
  ungroup() %>%
  group_by(parameter, model) %>%
  summarize(mean = round(mean(samples), 3),
            median = round(quantile(samples, .5), 3),
            x.mean = 1,
            y.mean = max(ifelse(str_detect(parameter, "beta"), 725, 
                                ifelse(str_detect(parameter, "delta"), 1950, 650))),
            x.median = 1,
            y.median = max(ifelse(str_detect(parameter, "beta"), 660, 
                                  ifelse(str_detect(parameter, "delta"), 1825, 575))))

corr.means <- corr.sims %>%
  ungroup() %>%
  group_by(parameter, model) %>%
  summarize(mean = round(mean(samples), 3),
            sd = round(sd(samples), 3),
            x.mean = -1,
            y.mean = max(ifelse(str_detect(parameter, "delta"), 1825, 625)),
            x.sd = -1,
            y.sd = max(ifelse(str_detect(parameter, "delta"), 1675, 575)))

sd.white.plot <- sd.sims %>%
  ggplot(aes(x = samples)) +
  geom_histogram(bins = 30) +
  facet_grid(parameter ~ model,
             scales = "free",
             labeller = label_parsed) +
  xlim(c(0, 2)) +
  theme_bw() +
  theme(text = element_text(family =  "LM Roman 10"))  +
  geom_text(data = sd.means, aes(x = x.mean, y = y.mean, label = paste0("Mean = ", mean)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  geom_text(data = sd.means, aes(x = x.median, y = y.median, label = paste0("Med. = ", median)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  xlab("SD") +
  ylab("Samples")


corr.white.plot <- corr.sims %>%
  ggplot(aes(x = samples)) + 
  geom_histogram(bins = 40) + 
  facet_grid(parameter ~ model,
             scales = "free",
             labeller = label_parsed) +
  xlim(c(-1, 1)) +
  theme_bw() +
  theme(text = element_text(family = "LM Roman 10"))  +
  geom_text(data = corr.means, aes(x = x.mean, y = y.mean, label = paste0("Mean = ", mean)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  geom_text(data = corr.means, aes(x = x.sd, y = y.sd, label = paste0("SD = ", sd)),
            size = 6 * .36, family = "LM Roman 10",
            hjust = 0) +
  xlab("Correlation") +
  ylab("Samples")

# ggsave(here("Figures", "sd_white_plot.pdf"), plot = sd.white.plot,
#        height = 5, width = 5, units = "in", device = "pdf",
#        dpi = 300)
# ggsave(here("Figures", "corr_white_plot.pdf"), plot = corr.white.plot,
#        height = 5, width = 5, units = "in", device = "pdf",
#        dpi = 300)



reglambda.sims <- make_mean_var_sims(rstan::extract(re3.reglambda, pars = c("theta0", "Sigma")))
lkj.sims <- make_mean_var_sims(rstan::extract(re3.lkj.stan, pars = c("theta0", "Sigma")))
iw.sims <- make_mean_var_sims(rstan::extract(re3.iw.stan, pars = c("theta0", "Sigma")))
cts0_from_one_posterior(reglambda.sims[1,])
reglambda.cts0.summary <- as.data.frame(cts_summary_from_gamma(reglambda.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS =  V1,
         Mean = V2,
         SD =   V3)
lkj.cts0.summary <- as.data.frame(cts_summary_from_gamma(lkj.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS = V1,
         Mean = V2,
         SD = V3)
iw.cts0.summary <- as.data.frame(cts_summary_from_gamma(iw.sims, M = 500)) %>%
  type_convert() %>%
  rename(CTS = V1,
         Mean = V2,
         SD = V3)

reglambda.cts0 <- reglambda.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
lkj.cts0 <- lkj.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
iw.cts0 <- iw.cts0.summary %>%
  transmute(CTS = CTS,
            Posterior = paste(round(Mean, 2), 
                              " (", round(`2.5%`, 2), 
                              ", ", round(`97.5%`, 2), ")", sep = ""))
all.cts0 <- cbind.data.frame(iw.cts0, lkj.cts0$Posterior, reglambda.cts0$Posterior)
names(all.cts0) <- c("CTS", "IW", "LKJ", "RHS-CS")

# print(xtable(all.cts0, type = "latex"), 
#       file = here("Results", "white-race-results.tex"), include.rownames = FALSE)

