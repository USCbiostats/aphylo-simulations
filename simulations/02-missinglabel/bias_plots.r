rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)

bias <- readRDS("simulations/02-missing/bias.rds")

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
sizelvls <- levels(bias$size_tag)
for (i in 1:4) {
  # Clearing plot space and creating the pdf
  
  pdf(sprintf("simulations/02-missing/bias_plots_tree-size=%s.pdf", switch (i,
    `1` = "small",
    `2` = "mid-small",
    `3` = "mid-large",
    `4` = "large"
  )), width = 6, height = 7)
  
  # Nobservations in this group
  nobs <- bias %>% dplyr::filter(as.integer(size_tag) == i) %>%
    subset(select=index) %>% unique %>% nrow
  
  
  # Creating the plot
  p <- bias %>% dplyr::filter(as.integer(size_tag) == i) %>%
    # Creating the boxplot
    ggplot(aes(parameter, bias)) +
    geom_violin(aes(fill = Prior)) +
    theme_minimal(base_family = "serif") +
    # Adding an horizontal line, and spliting my % of missings
    geom_hline(yintercept = 0, lty=2) +
    facet_grid(miss_tag ~ .) +
    geom_hline(yintercept = 2/20 - 2/40, lty = 3, col="red") +
    ylim(-.2,.2) + ylab("Bias") + xlab("") + 
    # labs(
    #   title = "Empirical Bias",
    #   subtitle = paste("Includes trees of size",levels(bias$size_tag)[i])
    # ) +
    scale_x_discrete(
      labels=c(
        mu0  = expression(mu[paste(0,1)]),
        mu1  = expression(mu[10]),
        psi0 = expression(psi[paste(0,1)]),
        psi1 = expression(psi[10]),
        eta0 = expression(eta[0]),
        eta1 = expression(eta[0]),
        Pi   = expression(pi)
      ),
      limits = c("mu0", "mu1", "psi0", "psi1", "Pi")
    ) +
    scale_color_grey(aesthetics = "fill")
  
  # Adding a title
  message(
    sprintf("Bias distribution for trees of size %s", levels(bias$size_tag)[i]),
    subtitle = sprintf("# of observations: %i", nobs)
  )
  
  # Printing it on screen (need to do that explicitly on a loop)
  print(p)
  message("Level ", levels(bias$size_tag)[i], " done.")
  
  dev.off()
  
}
