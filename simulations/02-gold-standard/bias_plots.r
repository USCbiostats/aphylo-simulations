rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)

load("simulations/02-gold-standard/bias.rda")

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
  
  pdf(sprintf("simulations/02-gold-standard/bias_trees_of_size_%s.pdf", switch (i,
    `1` = "small",
    `2` = "mid-small",
    `3` = "mid-large",
    `4` = "large"
  )))
  
  # Nobservations in this group
  nobs <- bias %>% dplyr::filter(as.integer(size_tag) == i) %>%
    subset(select=index) %>% unique %>% nrow
  
  
  # Creating the plot
  p <- bias %>% dplyr::filter(as.integer(size_tag) == i) %>%
    # Creating the boxplot
    ggplot(aes(parameter, bias)) + geom_boxplot(aes(colour = Prior)) +
    
    # Adding an horizontal line, and spliting my % of missings
    geom_hline(yintercept = 0, lty=2) + facet_grid(miss_tag ~ .) +
    ylim(-.2,.2) + ylab("Bias") + xlab("")
  
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
