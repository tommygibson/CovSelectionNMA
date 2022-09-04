###### Summarize results for Dong 2013

library(rstan)
library(tidyverse)
library(loo)
library(here)

dong_results <- read_rds(here("Results", "dong-models.rds"))
pagliaro_results <- read_rds(here("Results", "pagliaro-models.rds"))

dong_results[[2]]
round(summary(dong_results[[1]][[4]], pars = c("SD", "sigma_beta", "CORR"))$summary, 4)
