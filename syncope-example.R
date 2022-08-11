##### 3RE model with covariance selection for syncope data

library(rstan)
library(tidyverse)
library(here)

syncope <- read_csv(here("syncope.cleaned.csv")) %>%
  group_by(Variable) %>%
  filter(!is.na(n_i0),
         n() > 4)

male <- syncope %>%
  filter(Variable == "Male Gender")

re3.dat.stan <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                     a = 0, b = 2, c = 0, d = 2, f = 0, g = 2)
re3.dat.reghorse <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                         a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                         scale_global_lambda = 1, scale_global_gamma = 2,
                         nu_global_lambda = 1, nu_global_gamma = 3,
                         nu_local_lambda = 1, nu_local_gamma = 1,
                         slab_scale_lambda = 2, slab_df_lambda = 1,
                         slab_scale_gamma = 1, slab_df_gamma = 4)