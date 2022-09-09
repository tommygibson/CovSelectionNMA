##### Synthetic data example for 3RE models

library(rstan)
library(MASS)
library(here)


re3.std.model <- stan_model(here("Models", "3RE.stan"))
re3.iw.model <- stan_model(here("Models", "3RE-IW.stan"))
re3.lkj.model <- stan_model(here("Models", "3RE-LKJ.stan"))
re3.reglambda.model <- stan_model(here("Models", "3RE-regularized-lambda-cholesky.stan"))

# either 3, 5, or 7 studies (to mimic syncope data)
# S <- c(3, 5, 7)
set.seed(1234)
S <- 5

# keep prob of event for \RFbar group constant at 0.05
pi0 <- 0.05
# prob of event for \RF group varies between 0.075 and 0.2
# evenly spaced
# pi1 <- list()
# for(i in 1:length(S)){
#   pi1[[i]] <- seq(0.075, 0.2, length.out = S[i])
# }
pi1 <- seq(0.075, 0.3, length.out = S)

# prob of \RF is constant at 0.25
psi <- 0.25

# number of subjects per study will be uniform(250, 2500) as in 3RE paper


N <- sample(seq(500, 2500, 500), S, replace = FALSE)
pi11 <- pi1 * psi
pi10 <- (1 - pi1) * psi
pi01 <- pi0 * (1 - psi)
pi00 <- (1 - pi0) * (1 - psi)

probs <- cbind(pi11, pi10, pi01, pi00)
tabs <- matrix(nrow = S, ncol = 4)

for(i in 1:S){
  tabs[i,] <- rmultinom(1, N[i], probs[i,])
}


y <- tabs[,c(1, 3)]
n <- cbind(rowSums(tabs[, c(1, 2)]),
           rowSums(tabs[, c(3, 4)]))

nu.hat <- log((n[,1] / N) / (1 - n[,1] / N))
beta.hat <- apply(log((y/n) / (1 - y/n)), 1, mean)
delta.hat <- log((y[,1] * (n[,2] - y[,2])) / (y[,2] * (n[,1] - y[,1])))

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

re3.dat.reglambda <- list(S = S, y1 = y[,1], y0 = y[,2], n1 = n[,1], n0 = n[,2],
                          a = 0, b = 2, c = 0, d = 2, f = 0, g = 2,
                          scale_global_lambda = scale_global_lambda, scale_gamma = 2,
                          nu_global_lambda = 3,  nu_local_lambda = 1,
                          slab_scale_lambda = 1, slab_df_lambda = 4)

re3.stan <- sampling(re3.std.model, data = re3.dat.stan,
                     chains = 4, iter = 4000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     seed = 714)
re3.iw.stan <- sampling(re3.iw.model, data = re3.dat.stan,
                        chains = 4, iter = 4000,
                        pars = c("theta0", "SD", "CORR", "Sigma"),
                        seed = 714)

re3.lkj.stan <- sampling(re3.lkj.model, data = re3.dat.stan,
                         chains = 4, iter = 4000, # iter = 10000, warmup = 6000,
                         pars = c("theta0", "SD", "CORR", "Sigma"),
                         control = list(adapt_delta = .99, max_treedepth = 15),
                         cores = 4,
                         seed = 714)

re3.reglambda <- sampling(re3.reglambda.model, data = re3.dat.reglambda,
                          chains = 4, iter = 4000, warmup = 2000,
                          pars = c("theta0", "SD", "CORR", "Sigma"),
                          control = list(adapt_delta = 0.99, max_treedepth = 15),
                          cores = 4,
                          seed = 714)
round(summary(re3.stan, pars = c("beta0", "delta0", "nu0", "sigma_beta", "sigma_delta", "sigma_nu", "lp__"))$summary, 4)
round(summary(re3.iw.stan, pars = c("theta0", "SD", "CORR"))$summary, 4)
round(summary(re3.lkj.stan, pars = c("theta0", "SD", "CORR"))$summary, 4)
# round(summary(re3.reghorse)$summary, 4)
round(summary(re3.reglambda, pars = c("theta0", "SD", "CORR"))$summary, 4)

par(mfrow = c(3, 3))
for(i in 1:3){
  for(j in 1:3){
    hist(rstan::extract(re3.iw.stan, pars = "CORR")$CORR[,i,j], xlim = c(-1, 1), breaks = 30)
  }
}
for(i in 1:3){
  for(j in 1:3){
    hist(rstan::extract(re3.lkj.stan, pars = "CORR")$CORR[,i,j], xlim = c(-1, 1), breaks = 30)
  }
}
for(i in 1:3){
  for(j in 1:3){
    hist(rstan::extract(re3.reglambda, pars = "CORR")$CORR[,i,j], xlim = c(-1, 1), breaks = 30)
  }
}

sd.lim.max <- c(2, 2, 1.5)

sd.sims <- cbind.data.frame(c(as.vector(do.call(cbind, rstan::extract(re3.stan, pars = c("sigma_beta", "sigma_delta", "sigma_nu")))),
                              as.vector(rstan::extract(re3.iw.stan, pars = "SD")$SD),
                              as.vector(rstan::extract(re3.lkj.stan, pars = "SD")$SD),
                              as.vector(rstan::extract(re3.reglambda, pars = "SD")$SD)),
                            rep(rep(c("sigma_beta", "sigma_delta", "sigma_nu"), each = 8000), 4),
                            rep(c("std", "IW", "LKJ", "RHS"), each = 24000))
names(sd.sims) <- c("samples", "parameter", "model")

# extract only unique correlation values rather than whole 3x3 matrix
# only need rho_beta_delta, rho_beta_nu, rho_delta_nu, which are in the 2nd, 3rd, and 6th spots

corr.sims <- cbind.data.frame(c(as.vector(matrix(rstan::extract(re3.iw.stan, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)]),
                                as.vector(matrix(rstan::extract(re3.lkj.stan, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)]),
                                as.vector(matrix(rstan::extract(re3.reglambda, pars = "CORR")$CORR, 8000, 9)[, c(2, 3, 6)])),
                              rep(rep(c("rho_beta_delta", "rho_beta_nu", "rho_delta_nu"), each = 8000), 3),
                              rep(c("IW", "LKJ", "RHS"), each = 24000))
names(corr.sims) <- c("samples", "parameter", "model")

sd.sims %>%
  ggplot(aes(x = samples)) +
  geom_histogram(bins = 40) +
  facet_grid(parameter ~ model,
             scales = "free") +
  xlim(c(0, 1.5)) +
  theme_bw()

corr.sims %>%
  ggplot(aes(x = samples)) + 
  geom_histogram(bins = 40) + 
  facet_grid(parameter ~ model,
             scales = "free") +
  xlim(c(-1, 1)) +
  theme_bw()

for(i in 1:3){
  hist(rstan::extract(re3.iw.stan, pars = "SD")$SD[,i], breaks = 40, xlim = c(0, sd.lim.max[i]))
}
for(i in 1:3){
  hist(rstan::extract(re3.lkj.stan, pars = "SD")$SD[,i], breaks = 20, xlim = c(0, sd.lim.max[i]))
}
for(i in 1:3){
  hist(rstan::extract(re3.reghorse, pars = "SD")$SD[,i], breaks = 20, xlim = c(0, sd.lim.max[i]))
}


