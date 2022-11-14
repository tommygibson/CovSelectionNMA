#### testing different values for the conditional shrinkage prior

var.1 <- seq(0.01, 0.75, 0.02)
var.2 <- seq(0.01, 0.75, 0.02)

var.1.2 <- expand.grid(var.1, var.2)


shrinkage.prior <- function(x, a0){
  var1 <- x[1]
  var2 <- x[2]
  
  shrunk_var <- a0 * (1 / (1 / var1 + 1 / var2))
  return(shrunk_var)
}

var.1.2$shrunk_var1 <- apply(var.1.2[,c(1:2)], 1, shrinkage.prior, a0 = 1)
var.1.2$shrunk_var2 <- apply(var.1.2[,c(1:2)], 1, shrinkage.prior, a0 = 4)
var.1.2$shrunk_var3 <- apply(var.1.2[,c(1:2)], 1, shrinkage.prior, a0 = 9)
var.1.2$shrunk_var4 <- apply(var.1.2[,c(1:2)], 1, shrinkage.prior, a0 = 16)

var.1.2 %>%
  ggplot(aes(x = Var1, y = Var2, z = shrunk_var2)) +
  metR::geom_contour_fill() +
  geom_contour2(aes(label = stat(level)))
