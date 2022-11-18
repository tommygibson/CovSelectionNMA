#### testing different values for the conditional shrinkage prior

library(tidyverse)
library(metR)
library(here)

sd.1 <- seq(0.005, 1, 0.02)
sd.2 <- seq(0.005, 1, 0.02)

sd.1.2 <- expand.grid(sd.1, sd.2)


shrinkage.prior <- function(x, a0){
  sd1 <- x[1]
  sd2 <- x[2]
  
  shrunk_sd <- sqrt(a0 ^ 2 * (1 / (1 / (sd1 ^ 2) + 1 / (sd2 ^ 2))))
  return(shrunk_sd)
}

sd.1.2$shrunk_sd1 <- apply(sd.1.2[,c(1:2)], 1, shrinkage.prior, a0 = 1)
sd.1.2$shrunk_sd2 <- apply(sd.1.2[,c(1:2)], 1, shrinkage.prior, a0 = 2)
sd.1.2$shrunk_sd3 <- apply(sd.1.2[,c(1:2)], 1, shrinkage.prior, a0 = 3)
sd.1.2$shrunk_sd4 <- apply(sd.1.2[,c(1:2)], 1, shrinkage.prior, a0 = 4)

(cond.shrink.plot <- sd.1.2 %>%
  ggplot(aes(x = Var1, y = Var2, z = shrunk_sd4)) +
  #metR::geom_contour_fill() +
  #geom_contour() +
  geom_contour2(aes(label = stat(level)), breaks = c(0.2, 0.5, 0.75, 1, 1.5, 2)) +
  scale_color_viridis_c() +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  theme_bw() +
  labs(y = expression(omega[p]),
       x = expression(omega[q])) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)))
ggsave(here("Figures", "cond_shrink_plot.pdf"), plot = cond.shrink.plot,
       height = 4, width = 4, units = "in", device = "pdf", dpi = 300)
