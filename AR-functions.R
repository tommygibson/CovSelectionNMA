##### Functions for calculating ARs, LORs, etc from mean/variance
##### similar to calculating CTS_0 from mean/variance estimates in 3RE model

# function to calculate AR for each of T treatments from one posterior sample
# data structure (mean.sd) is (mu1, mu2,..., muT, sd1, ..., sdT)


get_AR <- function(mean.sd, M = 1000){
  ntrt <- length(mean.sd) / 2
  mu <- mean.sd[1:ntrt]
  sigma <- mean.sd[(ntrt + 1):(2 * ntrt)]
  
  pred <- mapply(function(x, y){1 / (1 + exp(-rnorm(M, mean = x, sd = y)))},
                 mu, sigma)
  
  return(apply(pred, 2, mean))
}
get_AR_SV <- function(mean.sd.beta, M = 1000){
  ntrt <- floor((length(mean.sd.beta) - 1) / 2)
  mu <- mean.sd.beta[1:ntrt]
  sigma <- mean.sd.beta[(ntrt + 1):(2 * ntrt)]
  sigma.beta <- mean.sd.beta[length(mean.sd.beta)]
  tot_sd <- sqrt(sigma ^ 2 + sigma.beta ^ 2)
  pred <- mapply(function(x, y){1 / (1 + exp(-rnorm(M, mean = x, sd = y)))},
                 mu, tot_sd)
  return(apply(pred, 2, mean))
}

get_full_AR <- function(mean.sd.mat, M = 1000){
  t(apply(X = mean.sd.mat, MARGIN = 1, FUN = get_AR, M = M))
}

get_full_AR_SV <- function(mean.sd.beta.mat, M = 1000){
  t(apply(X = mean.sd.beta.mat, MARGIN = 1, FUN = get_AR_SV, M = M))
}

simple_summary <- function(x){
  return(c(mean(x), sd(x), quantile(x, c(.025, .975))))
}
summarize_full_AR <- function(mean.sd.mat, M = 1000){
  
  full_AR <- get_full_AR(mean.sd.mat, M = M)
  
  if(is.null(dim(full_AR))){
    return(simple_summary(full_AR))
  } 
  
  else {
    return(t(apply(full_AR, 2, simple_summary)))
  }
}

summarize_full_AR_SV <- function(mean.sd.beta.mat, M = 1000){
  
  full_AR <- get_full_AR_SV(mean.sd.beta.mat, M = M)
  
  if(is.null(dim(full_AR))){
    return(simple_summary(full_AR))
  } 
  
  else {
    return(t(apply(full_AR, 2, simple_summary)))
  }
}

