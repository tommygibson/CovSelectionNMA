###### testing Copas-3RE

library(tidyverse)
library(R2jags)
library(here)
library(MASS)

# initial number of studies and subjects per study

set.seed(63)
S.init <- 100
N.init <- floor(runif(S.init, 100, 500))

# hyperparameters
beta0 <- log(0.3 / 0.7)
nu0 <- log(.15 / .85)
delta0 <- 3
sigma.beta <- 0.4
sigma.nu <- 0.4
sigma.delta <- 0.5
rho <- 0.7
gamma0 <- -0.3
gamma1 <- 0.06

mu <- c(beta0, delta0, nu0, 0)
Sigma <- rbind(c(sigma.beta^2, 0, 0, 0),
               c(0, sigma.delta^2, 0, rho * sigma.delta),
               c(0, 0, sigma.nu, 0),
               c(0, rho * sigma.delta, 0, 1))

RE <- mvrnorm(n = S.init, mu = mu, Sigma = Sigma)

tables <- as.data.frame(matrix(0, nrow = S.init, ncol = 5))
names(tables) <- c("TP", "FP", "FN", "TN", "sqrtESS")
ESS <- vector(length = S.init)

for(i in 1:S.init){
  beta <- RE[i, 1]
  delta <- RE[i, 2]
  nu <- RE[i, 3]
  epsilon <- RE[i, 4]
  pi1 <- 1 / (1 + exp(-(beta + delta / 2)))
  pi0 <- 1 / (1 + exp(-(beta - delta / 2)))
  psi <- 1 / (1 + exp(-nu))
  
  pi11 <- pi1 * psi
  pi10 <- (1 - pi1) * psi
  pi01 <- pi0 * (1 - psi)
  pi00 <- (1 - pi0) * (1 - psi)
  
  probs <- c(pi11, pi10, pi01, pi00)
  
  tables[i, 1:4] <- rmultinom(1, N.init[i], probs)
  
  n1 <- tables[i, 1] + tables[i, 3]  # events
  n2 <- tables[i, 2] + tables[i, 4]  # non-events
  
  tables[i, 5] <- 1 / sqrt(1 / n1 + 1 / n2)
  
}

dat.full <- cbind.data.frame(RE, tables, gamma0 + gamma1 * tables[,5] + RE[,4],
                             log(((tables[,1] + .5) * (tables[,4] + .5) / ((tables[,2] + .5) * (tables[,3] + .5)))),
                             sqrt(rowSums(1 / (tables[,1:4] + .5))))

names(dat.full)[c(1:4, 10:12)] <- c("beta.i", "delta.i", "nu.i", "epsilon.i", "z.i", "lnORhat", "SElnORhat")

dat.full %>%
  ggplot(aes(x = lnORhat, y = 1 / sqrtESS, color = (z.i > 0))) +
    geom_point()

dat.select <- dat.full %>%
  filter(z.i > 0) 

S <- nrow(dat.select)
y <- dat.select %>%
  dplyr::select(TP, FP)
n <- dat.select %>%
  transmute(n1 = TP + FN,
            n2 = FN + TN)
ESS <- dat.select$sqrtESS ^ 2
n.tot = rowSums(n)
a <- -1
b <- 0.1
c <- 0
d <- 0.1
e <- -1
f <- 0.1
B.beta <- B.delta <- B.nu <- 2

dat.copas3re <- list(S = S, y = y, n = n, n.tot = n.tot, ESS = ESS, 
                a = a, b = b, c = c, d = d, e = e, f = f,
                B.beta = B.beta, B.delta = B.delta, B.nu = B.nu)
dat.3re <- list(S = S, y = y, n = n, n.tot = n.tot, 
                a = a, b = b, c = c, d = d, e = e, f = f,
                B.beta = B.beta, B.delta = B.delta, B.nu = B.nu)

init.gen.copas3re <- function(){
  list(beta = runif(S, -1, 1),
       delta = runif(S, -1, 1),
       nu = runif(S, -1, 1),
       beta0 = runif(1, -1, 1),
       delta0 = runif(1, -1, 1),
       nu0 = runif(1, -1, 1),
       rho.delta.z = runif(1, -0.5, 0.5),
       sigma.beta = runif(1, 0.2, 1),
       sigma.delta = runif(1, 0.2, 1),
       sigma.nu = runif(1, 0.2, 1),
       z = runif(S, 0.1, 1),
       gamma0 = runif(1, -1, 1),
       gamma1 = runif(1, 0, 0.1)
  )
}
init.gen.copas3re.select <- function(){
  list(beta = runif(S, -1, 1),
       delta = runif(S, -1, 1),
       nu = runif(S, -1, 1),
       beta0 = runif(1, -1, 1),
       delta0 = runif(1, -1, 1),
       nu0 = runif(1, -1, 1),
       rho.spike = rbinom(1, 1, 0.5),
       rho.slab = runif(1, -0.5, 0.5),
       sigma.beta = runif(1, 0.2, 1),
       sigma.delta = runif(1, 0.2, 1),
       sigma.nu = runif(1, 0.2, 1),
       z = runif(S, 0.1, 1),
       gamma0 = runif(1, -1, 1),
       gamma1 = runif(1, 0, 0.1)
  )
}

params.copas3re <- c("delta0", "rho.delta.z", "sigma.delta")

fit.copas3re <- do.call(jags.parallel,
                        list(data = names(dat.copas3re), init.gen.copas3re, params.copas3re,
                             model.file = here("Models", "Copas-3RE.txt"),
                             n.chains = 4, n.iter = 100000, n.burnin = 40000, n.thin = 8, DIC = FALSE))
fit.copas3re.select <- do.call(jags.parallel,
                        list(data = names(dat.copas3re), init.gen.copas3re.select, params.copas3re,
                             model.file = here("Models", "Copas-3RE-select.txt"),
                             n.chains = 4, n.iter = 100000, n.burnin = 40000, n.thin = 8, DIC = FALSE))

fit.copas3re$BUGSoutput$summary
fit.copas3re.select$BUGSoutput$summary
hist(fit.copas3re$BUGSoutput$sims.list$rho.delta.z)
hist(fit.copas3re.select$BUGSoutput$sims.list$rho.delta.z)




################# Testing regular copas vs copas-select
set.seed(77)
S.full <- 50
theta0 <- 0.5
tau <- 0.1
gamma0 <- -0.5
gamma1 <- 0.3
rho <- 0.9
s.full <- runif(S.full, 0.05, 0.8)
theta.full <- rnorm(S.full, theta0, tau)
y.z <- as.data.frame(matrix(nrow = S.full, ncol = 3))
names(y.z) <- c("y", "z", "s")
for(i in 1:S.full){
  y.z[i, 1:2] <- mvrnorm(1, mu = c(theta.full[i], gamma0 + gamma1 * (1 / s.full[i])), 
                 Sigma = matrix(c(s.full[i]^2, rho * s.full[i], 
                                rho * s.full[i], 1), nrow = 2, byrow = T))
}
y.z$s <- s.full

y.z.select <- y.z %>%
  filter(z > 0)

S <- nrow(y.z.select)
y <- y.z.select$y
s <- y.z.select$s

# S <- nrow(y.z)
# y <- y.z$y
# s <- y.z$s

dat.copas <- list(y = y,
                  s = s,
                  S = S)
dat.copas.select <- list(y = y,
                         s = s,
                         S = S,
                         p = 0.5)

init.gen.copas <- function(){
  list(theta0 = runif(1, -1, 1),
       tau = runif(1, 0.2, 1),
       rho = runif(1, -0.5, 0.5)
  )
}
init.gen.copas.select <- function(){
  list(theta0 = runif(1, -1, 1),
       tau = runif(1, 0.2, 1),
       rho.slab = runif(1, -0.5, 0.5),
       rho.spike = rbinom(1, 1, 0.5))
}
params.copas <- c("theta0", "tau", "rho", "gamma0", "gamma1")
params.copas.select <- c("theta0", "tau", "rho", "rho.spike", "rho.slab", "gamma0", "gamma1")

fit.copas <- jags(data = dat.copas, inits = init.gen.copas, parameters.to.save = params.copas,
                  model.file = here("Models", "Copas-std.txt"),
                  n.iter = 40000, n.burnin = 10000, n.chains = 4, n.thin = 4, DIC = FALSE)
fit.copas.select <- jags(data = dat.copas.select, inits = init.gen.copas.select, parameters.to.save = params.copas.select,
                         model.file = here("Models", "Copas-std-select.txt"),
                         n.iter = 40000, n.burnin = 10000, n.chains = 4, n.thin = 4, DIC = FALSE)


