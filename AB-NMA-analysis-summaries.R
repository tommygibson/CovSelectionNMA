###### Summarize results for Dong 2013

library(rstan)
library(tidyverse)
library(loo)
library(here)
library(xtable)
library(metadat)
library(gridExtra)
library(extrafont)

source(here("lay_out.R"))

source(here("AR-functions.R"))

# probably unnecessary to fully load the data
nma.dat.1 <- read_csv(here("Data", "dat.dong2013.csv")) 

len <- dim(nma.dat.1)[1]
ntrt <- length(unique(nma.dat.1$treatment))
nstudy <- max(nma.dat.1$id)

nma.dat.1$treatment <- factor(nma.dat.1$treatment,
                              levels = names(table(nma.dat.1$treatment))[order(table(nma.dat.1$treatment))])
dong_6 <- read_rds(here("Results" , "dong-models6-allarms.rds"))[[1]]

# Results from Dong 2013 analysis wihtout treatment arm w/ only 2 observations 

# dong_models_6arms <- read_rds(here("Results", "dong-models1-5-allarms-ordered.rds"))[[1]]
# dong_loos_6arms <- read_rds(here("Results", "dong-models1-5-allarms-ordered.rds"))[[2]]
# # add model 6, which ran separately
# 
# dong_models_6arms[[6]] <- read_rds(here("Results", "dong-models6-allarms.rds"))[[1]]
# dong_loos_6arms[[6]] <- read_rds(here("Results", "dong-models6-allarms.rds"))[[2]]
# 
# names(dong_models_6arms)[6] <- names(dong_loos_6arms)[6] <- "RHS-NS-SV"
# 
# round(summary(dong_models_6arms[[6]], pars = c("MU", "AR", "lp__"))$summary, 4)
# round(summary(dong_models_6arms[[1]], pars = c("SD"))$summary, 4)
# 
# waic(extract_log_lik(dong_models_6arms[[1]], parameter_name = "log_likelihood"))
# waic(extract_log_lik(dong_models_6arms[[6]], parameter_name = "log_likelihood"))

dong_models <- read_rds(here("Results", "dong-models-ordered.rds"))[[1]]
dong_loos <- read_rds(here("Results", "dong-models-ordered.rds"))[[2]]

#### Table of ELPDs

# elpds_6arms <- as.data.frame(matrix(nrow = length(dong_loos_6arms), ncol = 3))
elpds <- as.data.frame(matrix(nrow = length(dong_loos), ncol = 3))
# elpds_6arms[,1] <- c("IW", "LKJ", "RHS-CS", "RHS-NS", "RHS-CS-SV", "RHS-NS-SV")
elpds[, 1] <- c("IW", "LKJ", "RHS-CS", "RHS", "RHS-SV")
for(i in 1:length(dong_loos)){
  elpds[i, 2:3] <- dong_loos[[i]]$estimates[c(1, 3), 1]
}

names(elpds) <- c("Model", "ELPD", "LOOIC")

print(xtable(elpds, type = "latex", digits = 2),
      file = here("Results", "SIM-elpds.tex"), include.rownames = FALSE)


###### Plot of dong data, showing that it's basically just study-level variation

dat.dong2013 %>%
  mutate(logodds = log((death + 0.5) / (randomized + 1)) - log((randomized - death + 0.5) / (randomized + 1))) %>%
  ggplot(aes(x = id, y = logodds, color = treatment)) + 
  geom_point() +
  theme_bw() +
  xlab("Study") +
  ylab("log-odds of event") +
  scale_color_discrete(name = "Treatment")


###### Table of AR posteriors for each model

#nsims <- dim(rstan::extract(dong_models_6arms[[1]], pars = "SD")$SD)[1]

#corr.indices <- c(2:6, 9:12, 16:18, 23:24, 30)
treatment.order <- levels(nma.dat.1$treatment)

param.summaries.list <- list()
for(i in 1:length(dong_models)){

  # one row per parameter (SD, AR, MU)
  param.summaries.list[[i]] <- as.data.frame(matrix(nrow = 6 + 6 + 6, ncol = 7))
  names(param.summaries.list[[i]]) <- c("Mean", "SD", "lower.ci", "upper.ci", "Parameter", "Treatment", "Model")
  # standard deviations
  param.summaries.list[[i]][1:6, 1:4] <- summary(dong_models[[i]], pars = "SD")$summary[, c(1, 3, 4, 8)]
  # absolute risks
  param.summaries.list[[i]][7:12, 1:4] <- summary(dong_models[[i]], pars = "AR")$summary[1:6, c(1, 3, 4, 8)]
  # treatment mean log-odds
  param.summaries.list[[i]][13:18, 1:4] <- summary(dong_models[[i]], pars = "MU")$summary[1:6, c(1, 3, 4, 8)]
  # column containing what parameter the mean/sd/CI correspond to
  param.summaries.list[[i]][, 5] <- c(paste0("SD", 1:6),
                                      paste0("AR", 1:6),
                                      paste0("MU", 1:6))
  # column containing treatments
  param.summaries.list[[i]][, 6] <- c(levels(nma.dat.1$treatment),
                                      rep(levels(nma.dat.1$treatment), 2))

  # column containing what model the estimates came from
  param.summaries.list[[i]][, 7] <- names(dong_models)[i]
}

# for(i in 1:length(dong_models_6arms)){
#   
#   # one row per parameter (SD, Correlations, AR, MU)
#   # param.summaries.list[[i]] <- as.data.frame(matrix(nrow = 6 + (6 * 5 / 2) + 6 + 6, ncol = 7))
#   param.summaries.list[[i]] <- as.data.frame(matrix(nrow = 6 + 6 + 6, ncol = 7))
#   names(param.summaries.list[[i]]) <- c("Mean", "SD", "lower.ci", "upper.ci", "Parameter", "Treatment", "Model")
#   # standard deviations
#   param.summaries.list[[i]][1:6, 1:4] <- summary(dong_models_6arms[[i]], pars = "SD")$summary[, c(1, 3, 4, 8)]
#   # correlations
#   #param.summaries.list[[i]][7:21, 1:4] <- summary(dong_models_6arms[[i]], pars = "CORR")$summary[corr.indices, c(1, 3, 4, 8)]
#   # absolute risks
#   param.summaries.list[[i]][22:27, 1:4] <- summary(dong_models_6arms[[i]], pars = "AR")$summary[1:6, c(1, 3, 4, 8)]
#   # treatment mean log-odds
#   param.summaries.list[[i]][28:33, 1:4] <- summary(dong_models_6arms[[i]], pars = "MU")$summary[1:6, c(1, 3, 4, 8)]
#   # column containing what parameter the mean/sd/CI correspond to
#   param.summaries.list[[i]][, 5] <- c(paste0("SD", 1:6), 
#                                       paste0("CORR", c(paste0(1, 2:6), paste0(2, 3:6), paste0(3, 4:6), paste0(4, 5:6), "56")),
#                                       paste0("AR", 1:6),
#                                       paste0("MU", 1:6))
#   # column containing treatments
#   param.summaries.list[[i]][, 6] <- c(levels(nma.dat.1$treatment), 
#                                       rep(NA, 15),
#                                       rep(levels(nma.dat.1$treatment), 2))
#                                       
#   # column containing what model the estimates came from
#   param.summaries.list[[i]][, 7] <- names(dong_models_6arms)[i]
# }

# merge each model's estimates
param.summaries <- do.call(rbind, param.summaries.list)

# plot will have 5 treatments on left side, and TIO-SMI on right
AR1 <- param.summaries %>%
  filter(str_detect(Parameter, "^AR"),
         !str_detect(Parameter, "1$")) %>%
  ggplot(aes(x = Mean, y = Model, color = Model, shape = Model)) + 
  geom_point() +
  geom_linerange(aes(xmin = lower.ci, xmax = upper.ci, y = Model, color = Model)) +
  scale_shape_manual(values = c(0, 2, 4, 15, 17)) +
  facet_wrap(Treatment ~., nrow = 5) +
  theme_bw() +
  theme(text = element_text("LM Roman 10"),
        #legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab(NULL) +
  xlab("Absolute Risk")

AR2 <- param.summaries %>%
  filter(str_detect(Parameter, "AR1")) %>%
  ggplot(aes(x = Mean, y = Model, color = Model, shape = Model)) + 
  geom_point() +
  geom_linerange(aes(xmin = lower.ci, xmax = upper.ci, y = Model, color = Model)) +
  scale_shape_manual(values = c(0, 2, 4, 15, 17)) +
  facet_wrap(Treatment ~., nrow = 1) +
  theme_bw() +
  theme(text = element_text("LM Roman 10"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.ticks.y = element_blank()) +
  ylab(NULL) +
  xlab("Absolute Risk")

pdf(file = here("Results", "dong-example-AR.pdf"),
    height = 6, width = 6)

AR1.2 <- grid.arrange(AR2, AR1, nrow = 1)

dev.off()

# ggsave("dong-example-AR.pdf", plot = AR1.2,
#        device = "pdf", height = 6, width = 5, units = "in")



##### Doing AR the right way (MC method)

iw.full <- do.call(cbind, rstan::extract(dong_models[[1]], pars = c("MU", "SD")))

ARs.iw <- summarize_full_AR(do.call(cbind, rstan::extract(dong_models[[1]], pars = c("MU", "SD"))))
ARs.lkj <- summarize_full_AR(do.call(cbind, rstan::extract(dong_models[[2]], pars = c("MU", "SD"))))
ARs.RHS.CS <- summarize_full_AR(do.call(cbind, rstan::extract(dong_models[[3]], pars = c("MU", "SD"))))
ARs.RHS <- summarize_full_AR(do.call(cbind, rstan::extract(dong_models[[4]], pars = c("MU", "SD"))))
ARs.RHS.SV <- summarize_full_AR_SV(do.call(cbind, rstan::extract(dong_models[[5]], pars = c("MU", "SD", "sigma_beta"))))

ARs.all <- cbind.data.frame(rbind(ARs.iw, ARs.lkj, ARs.RHS.CS,
                                  ARs.RHS, ARs.RHS.SV),
                            rep(levels(nma.dat.1$treatment), 5),
                            rep(c("IW", "LKJ", "RHS-CS", "RHS", "RHS-SV"), each = 6))

names(ARs.all) <- c("Mean", "SD", "lower.ci", "upper.ci", "Treatment", "Model")
rownames(ARs.all) <- NULL

ARs.all <- ARs.all %>% type_convert()

AR1 <- ARs.all %>%
  filter(Treatment != "TIO-SMI") %>%
  ggplot(aes(x = Mean, y = Model, color = Model, shape = Model)) + 
  geom_point() +
  geom_linerange(aes(xmin = lower.ci, xmax = upper.ci, y = Model, color = Model)) +
  scale_shape_manual(values = c(0, 2, 4, 15, 17)) +
  facet_wrap(Treatment ~., nrow = 5) +
  theme_bw() +
  theme(text = element_text("LM Roman 10"),
        #legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab(NULL) +
  xlab("Absolute Risk")

AR2 <- ARs.all %>%
  filter(Treatment == "TIO-SMI") %>%
  ggplot(aes(x = Mean, y = Model, color = Model, shape = Model)) + 
  geom_point() +
  geom_linerange(aes(xmin = lower.ci, xmax = upper.ci, y = Model, color = Model)) +
  scale_shape_manual(values = c(0, 2, 4, 15, 17)) +
  facet_wrap(Treatment ~., nrow = 1) +
  theme_bw() +
  theme(text = element_text("LM Roman 10"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.ticks.y = element_blank()) +
  ylab(NULL) +
  xlab("Absolute Risk")

pdf(file = here("Results", "dong-example-AR-MC.pdf"),
    height = 6, width = 6)

AR1.2 <- grid.arrange(AR2, AR1, nrow = 1)

dev.off()
