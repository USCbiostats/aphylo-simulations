rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(magrittr)

bias <- readRDS("simulations/01-gold-standard/bias.rds")
conv <- readRDS("simulations/01-gold-standard/gelman.rds")

# Filtering according to convergence
bias <- bias[conv$mpsrf <= 1.1, ]

bias <- bias %>%
  as_tibble %>%
  filter(Missing < 1, PropOf0 > 0, PropOf0 < 1) %>%
  select(
    index, Prior, size_tag, miss_tag, PropLeafs_tag, matches("_bias$"),
    pscore, pscore_rand, pscore_worse) %>%
  gather(
    key = "parameter",
    value = "bias",
    which(grepl("(?<=bias)$", colnames(.), perl = TRUE))
    ) %>%
  mutate(parameter = gsub("_bias", "", parameter)) 
  

#+ plotting, echo=TRUE
# Plot -------------------------------------------------------------------------
graphics.off()

# Clearing plot space and creating the pdf

pdf("simulations/01-gold-standard/bias_plots_01.pdf", width = 6, height = 4)

# Nobservations in this group
nobs <- bias %>% 
  subset(select=index) %>% unique %>% nrow

# Creating the plot
parameter_labels <- c(
  mu_d0  = "mu[paste(d0,1)]",
  mu_d1  = "mu[d10]",
  mu_s0  = "mu[paste(s0,1)]",
  mu_s1  = "mu[s10]",
  psi0 = "psi[paste(0,1)]",
  psi1 = "psi[10]",
  eta0 = "eta[0]",
  eta1 = "eta[0]",
  Pi   = "pi"
)

bias$parameter <- parameter_labels[bias$parameter]

p <-
  bias %>% # dplyr::filter(as.integer(size_tag) == i) %>%
  # Creating the boxplot
  ggplot(aes(bias, y=size_tag)) +
  geom_density_ridges(aes(fill = Prior), color = "white") + 
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  # Adding an horizontal line, and spliting my % of missings
  facet_grid(Prior ~ parameter, labeller = labeller(parameter=label_parsed)) +
  coord_flip() +
  geom_vline(xintercept = 0, lty=2) + 
  # geom_vline(xintercept = 2/20 - 2/40, lty = 3, col="red") +
  xlim(-.5,.5) +
  xlab("Bias") +
  ylab("Number of Leaves") +
  scale_fill_grey()
    

# Printing it on screen (need to do that explicitly on a loop)
print(p)

dev.off()
 
# Prediction score -------------------------------------------------------------
library(tidyr)
bias <- readRDS("simulations/01-gold-standard/bias.rds")

# Filtering according to convergence
bias <- bias[conv$mpsrf <= 1.1, ]

bias <- bias %>%
  as_tibble %>%
  filter(Missing < 1, PropOf0 > 0, PropOf0 < 1) %>%
  mutate(
    pscore      = 1 - pscore/pscore_worse,
    pscore_rand = 1 - pscore_rand/pscore_worse,
  ) %>%
  select(
    index, Prior, size_tag, miss_tag, PropLeafs_tag,
    pscore, pscore_rand, auc) %>%
  gather(
    key = "Type",
    value = "Score",
    starts_with("pscore"), auc
  ) %>% 
  mutate(
    Type = if_else(
      Type == "pscore",
      "P. Score Obs.",
      ifelse(Type == "auc", "AUC", "P. Score Rand."))
  )

ggplot(bias, aes(y=Score, x=Type)) +
  geom_violin() +
  theme_minimal(base_family = "serif") +
  facet_grid(~size_tag) + 
  ylab("Prediction Score") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  xlab("") 

ggsave("simulations/01-gold-standard/prediction_scores-01.pdf", width=6, height=4)
  
