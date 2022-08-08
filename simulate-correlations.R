set.seed(711)
S <- 7
N <- 1000
beta0 <- nu0 <- log(.15 / .85)
delta0 <- 2
rho.beta.delta <- 0.7
rho.beta.nu <- 0
rho.delta.nu <- 0.7
sigma.beta <- 0.5
sigma.delta <- 1
sigma.nu <- 0.01

Lambda <- diag(c(.5, .5, .01))
Gamma <- diag(3)
Gamma[2,1] <- 0.8
Gamma[3,1] <- 0
Gamma[3,2] <- 0.8
hyper.var <- Lambda %*% Gamma %*% t(Gamma) %*% Lambda
cormat <- matrix(nrow = 3, ncol = 3)
for(i in 1:3){
  for(j in 1:3){
    cormat[i, j] <- hyper.var[i, j] / (sqrt(hyper.var[i, i]) * sqrt(hyper.var[j, j]))
  }
}

hyper.mean <- c(beta0, delta0, nu0)
# hyper.var <- matrix(c(sigma.beta ^ 2, rho.beta.delta * sigma.beta * sigma.delta, rho.beta.nu * sigma.beta * sigma.nu,
#                       rho.beta.delta * sigma.beta * sigma.delta, sigma.delta ^ 2, rho.delta.nu * sigma.delta * sigma.nu,
#                       rho.beta.nu * sigma.beta * sigma.nu, rho.delta.nu * sigma.delta * sigma.nu, sigma.nu ^ 2), 
#                     nrow = 3, byrow = TRUE)
cor.hat <- matrix(nrow = 1000, ncol = 3)
for(i in 1:1000){
  
  theta <- mvrnorm(S, hyper.mean, hyper.var)
  beta <- theta[,1]
  delta <- theta[,2]
  nu <- theta[,3]
  
  pi1 <- 1 / (1 + exp(-(beta + delta / 2)))
  pi0 <- 1 / (1 + exp(-(beta - delta / 2)))
  psi <- 1 / (1 + exp(-nu))
  
  pi11 <- pi1 * psi
  pi10 <- (1 - pi1) * psi
  pi01 <- pi0 * (1 - psi)
  pi00 <- (1 - pi0) * (1 - psi)
  
  probs <- cbind(pi11, pi10, pi01, pi00)
  
  tabs <- t(apply(probs, 1, function(x) {
    rmultinom(1, N, prob = x)
  }))
  
  y <- tabs[,c(1, 3)]
  n <- cbind(rowSums(tabs[, c(1, 2)]),
             rowSums(tabs[, c(3, 4)]))
  bigN <- rep(N, S)
  
  nu.hat <- log((n[,1] / bigN) / (1 - n[,1] / bigN))
  beta.hat <- apply(log((y/n) / (1 - y/n)), 1, mean)
  delta.hat <- log((tabs[,1] * tabs[,4]) / (tabs[,2] * tabs[,3]))
  cor.hat[i,1] <- cor(beta.hat, delta.hat)
  cor.hat[i,2] <- cor(beta.hat, nu.hat)
  cor.hat[i,3] <- cor(delta.hat, nu.hat)
}