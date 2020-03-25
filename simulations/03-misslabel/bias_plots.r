rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(magrittr)

source("simulations/00-bias-functions.r")
bias <- readRDS("simulations/03-misslabel/bias.rds")
conv <- readRDS("simulations/03-misslabel/gelman.rds")
# Filtering according to convergence
bias <- bias[conv$mpsrf <= 1.1, ]

bias$miss_tag <- interval_tags(
  bias$Missing, c(.1, .3, .5, .7, .9, 1))

bias <- bias %>%
  as_tibble %>%
  filter(Missing < 1, PropOf0 > 0, PropOf0 < 1) %>%
  select(index, Prior, size_tag, miss_tag, PropLeafs_tag, matches("_bias$")) %>%
  gather(
    key = "parameter",
    value = "bias",
    which(grepl("(?<=bias)$", colnames(.), perl = TRUE))
  ) %>%
  mutate(parameter = gsub("_bias", "", parameter)) 

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


bias$parameter <- factor(
  match(bias$parameter, names(parameter_labels)),
  levels = 1:length(parameter_labels),
  labels = parameter_labels)

#+ plotting, echo=TRUE
# Plot -------------------------------------------------------------------------
graphics.off()

p_wrong <-
  bias %>% dplyr::filter(Prior == "Wrong") %>%
  # Creating the boxplot
  ggplot(aes(bias, y = miss_tag)) +
  geom_density_ridges(color = "white") + 
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  # Adding an horizontal line, and spliting my % of missings
  facet_grid( ~ parameter, labeller = labeller(parameter=label_parsed)) +
  coord_flip() +
  geom_vline(xintercept = 0, lty=2) + 
  xlim(-.5,.5) +
  xlab("Bias") +
  ylab("Proportion of Missing Labels") +
  scale_fill_grey()

ggsave("simulations/03-misslabel/bias_plots_03_wrong_prior.pdf", width = 7, height = 4)

p_right <-
  bias %>% dplyr::filter(Prior == "Right") %>%
  # Creating the boxplot
  ggplot(aes(bias, y = miss_tag)) +
  geom_density_ridges(color = "white") + 
  theme_minimal(base_family = "serif") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  # Adding an horizontal line, and spliting my % of missings
  facet_grid( ~ parameter, labeller = labeller(parameter=label_parsed)) +
  coord_flip() +
  geom_vline(xintercept = 0, lty=2) + 
  xlim(-.5,.5) +
  xlab("Bias") +
  ylab("Proportion of Missing Labels") +
  scale_fill_grey()

ggsave("simulations/03-misslabel/bias_plots_03_right_prior.pdf", width = 7, height = 4)

# Prediction score -------------------------------------------------------------
library(tidyr)
bias <- readRDS("simulations/03-misslabel/bias.rds")

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

ggplot(bias, aes(Score, y = miss_tag)) +
  geom_density_ridges() +
  facet_grid(Prior ~ Type) +
  theme_minimal(base_family = "serif") +
  xlab("Prediction Score") +
  ylab("") +
  xlim(0, 1) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  coord_flip()

ggsave("simulations/03-misslabel/prediction_scores-03.pdf", width=6, height=6)
