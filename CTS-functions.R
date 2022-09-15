# functions to calculate CTS0's from mean and covariance posteriors
library(MASS)

# function to take list of samples from stan object and put it in the format
# for calculating CTSs
make_mean_var_sims <- function(sims.list){
  theta <- sims.list$theta0
  Sigma.list <- list()
  for(i in 1:dim(sims.list$Sigma)[1]){
    Sigma.list[[i]] <- as.vector(sims.list$Sigma[i,,])
  }
  Sigma <- do.call(rbind, Sigma.list)
  return(cbind(theta, Sigma))
}

# calculate CTS0(paramters) for one posterior sample
cts0_from_one_posterior <- function(mean.var, M = 1000){
  
  mu <- mean.var[1:3]
  Sigma <- matrix(mean.var[4:12], nrow = 3, byrow = FALSE)
  
  theta.new <- mvrnorm(M, mu = mu, Sigma = Sigma)
  beta.new <- theta.new[, 1]
  delta.new <- theta.new[, 2]
  nu.new <- theta.new[, 3]
  
  pi1.new <- 1 / (1 + exp(-(beta.new + delta.new / 2)))
  pi0.new <- 1 / (1 + exp(-(beta.new - delta.new / 2)))
  psi.new <- 1 / (1 + exp(-(nu.new)))
  
  sens.new <- (pi1.new * psi.new) / (pi1.new * psi.new + pi0.new * ( 1 - psi.new))
  spec.new <- ((1 - pi0.new) * (1 - psi.new)) / ((1 - pi0.new) * (1 - psi.new) + (1 - pi1.new) * psi.new)
  LRm.new <- (1 - sens.new) / spec.new
  LRp.new <- sens.new / (1 - spec.new)
  
  mc.est <- apply(cbind(LRm.new, LRp.new, 1 - pi0.new, pi1.new, sens.new, spec.new), MARGIN = 2, mean)
  return(mc.est)
}

# calculate CTS0(parameters) for all samples
make_cts0_from_all_posterior <- function(hyperparameter_matrix, M = 1000){
  t(apply(hyperparameter_matrix, MARGIN = 1, cts0_from_one_posterior, M = M))
}

# make a summary of each CTS with mean, sd, quantiles
simple_summary <- function(vec){
  return(c(mean(vec, na.rm=T), sd(vec, na.rm=T), quantile(vec, c(.025, .5, .975), na.rm=T)))
}

# put it all together, plus throwing trimming problematic samples
cts_summary_from_gamma <- function(hyperparameter_matrix, M){
  
  temp_all_cts0 <- make_cts0_from_all_posterior(hyperparameter_matrix, M = M)
  bad_LRm <- temp_all_cts0[,1] > 2.5
  bad_LRp <- temp_all_cts0[,2] > 30
  temp_all_cts0 <- temp_all_cts0[!(bad_LRp | bad_LRm),]
  summary_0 <- t(apply(temp_all_cts0, MARGIN = 2, simple_summary))
  
  rownames(summary_0) <- NULL
  
  return(cbind(rep(c("LRm", "LRp", "NPV", "PPV", "sens", "spec")),
               summary_0))
  
}
