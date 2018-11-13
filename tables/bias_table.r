rm(list = ls())

# Setting xtable options
options(xtable.booktabs = TRUE)
options(xtable.include.rownames = FALSE)
options(xtable.sanitize.text.function = function(x) x)

library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(xtable)

load("simulations/bias.rda")

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


# Bias per Prior, Size and Missing
bias_table <- bias %>%
  group_by(parameter, Prior, size_tag, miss_tag) %>%
  summarize(
    lower = quantile(bias, .025),
    med   = quantile(bias, .500),
    upper = quantile(bias, .975)
  ) %>% ungroup %>%
  transmute(
    parameter,
    Prior,
    size_tag = paste0("{", size_tag,"}"),
    miss_tag = paste0("{", miss_tag,"}"),
    range_and_med = sprintf("{[%.2f, %.2f, %.2f]}", lower, med, upper)
  ) %>%
  spread(parameter, range_and_med) %>%
  rename(
    Missing = miss_tag,
    Size    = size_tag
  ) 
  
fact <- bias_table$Prior

bias_table %>%
  select(-Prior) %>%
  split(., fact) %>%
  `attr<-`("subheadings", paste(names(.), "Prior")) %>%
  xtableList(
    caption = "Bias statistics per type of prior, size of the tree and proportion of missings. Each set of square brackets contains the 2.5, 50, and 97.5 quantiles of the distribution.",
    label   = "tab:bias-prior-size-missigness") %>%
  print %>%
  cat(file = "tables/bias_by_method_size_missingness.tex")
  

# Bias per Prior, and Missing
bias_table <- bias %>%
  group_by(parameter, Prior, miss_tag) %>%
  summarize(
    lower = quantile(bias, .025),
    med   = quantile(bias, .500),
    upper = quantile(bias, .975)
  ) %>% ungroup %>%
  transmute(
    parameter,
    Prior,
    miss_tag = paste0("{", miss_tag,"}"),
    range_and_med = sprintf("{[%.2f, %.2f, %.2f]}", lower, med, upper)
  ) %>%
  spread(parameter, range_and_med) %>%
  rename(
    Missing = miss_tag
  ) 

fact <- bias_table$Prior

bias_table %>%
  select(-Prior) %>%
  split(., fact) %>%
  `attr<-`("subheadings", paste(names(.), "Prior")) %>%
  xtableList(
    caption = "Bias statistics per type of prior, and proportion of missings. Each set of square brackets contains the 2.5, 50, and 97.5 quantiles of the distribution.",
    label   = "tab:bias-prior-missigness") %>%
  print %>%
  cat(file = "tables/bias_by_method_missingness.tex")

