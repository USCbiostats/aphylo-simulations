rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)

bias <- readRDS("simulations/01-gold-standard/bias.rds")

bias <- select(
  bias, index, Prior, size_tag, miss_tag, PropLeafs_tag,
  matches("_bias$")
  ) %>%
  as_tibble %>%
  gather(
    key = "parameter",
    value = "bias",
    which(grepl("(?<=bias)$", colnames(.), perl = TRUE))
    ) %>%
  mutate(parameter = gsub("_.+", "", parameter)) 
  

#+ plotting, echo=TRUE
# Plot -------------------------------------------------------------------------
graphics.off()

# Clearing plot space and creating the pdf

pdf("simulations/01-gold-standard/bias_plots_tree.pdf", width = 6, height = 7)

# Nobservations in this group
nobs <- bias %>% 
  subset(select=index) %>% unique %>% nrow

# Creating the plot
p <-
  bias %>% # dplyr::filter(as.integer(size_tag) == i) %>%
  # Creating the boxplot
  ggplot(aes(parameter, bias)) +
  geom_violin(aes(fill = Prior)) +
  theme_minimal(base_family = "serif") +
  # Adding an horizontal line, and spliting my % of missings
  geom_hline(yintercept = 0, lty=2) + 
  geom_hline(yintercept = 2/20 - 2/40, lty = 3, col="red") +
  facet_grid(size_tag ~ .) +
  ylim(-.2,.2) + ylab("Bias") + xlab("") + 
  scale_x_discrete(
    labels=c(
      mu0  = expression(mu[paste(0,1)]),
      mu1  = expression(mu[10]),
      psi0 = expression(psi[paste(0,1)]),
      psi1 = expression(psi[10]),
      eta0 = expression(eta[0]),
      eta1 = expression(eta[0]),
      Pi   = expression(pi)
    )
  ) +
    scale_color_grey(aesthetics = "fill")
    

# Printing it on screen (need to do that explicitly on a loop)
print(p)

dev.off()
 
