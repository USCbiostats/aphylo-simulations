rm(list = ls())

library(dplyr)
library(magrittr)
library(xtable)

# Setting xtable options
options(xtable.booktabs = TRUE)
options(xtable.include.rownames = FALSE)
options(xtable.sanitize.text.function = function(x) x)

load("simulations/03_uniform_drop_of_0/bias.rda")

# Cleaning latex errors
# See https://www.sharelatex.com/learn/Errors/Illegal_unit_of_measure_(pt_inserted)
bias <- bias %>%
  mutate(
    size_tag = paste0("{", size_tag, "}"),
    miss_tag = paste0("{", miss_tag, "}")
  )

# Coverage probability
summ <- bias %>% 
  group_by(Prior, size_tag, miss_tag) %>%
  summarize(
    psi0 = mean(psi0_covered95),
    psi1 = mean(psi1_covered95),
    mu0 = mean(mu0_covered95),
    m1 = mean(mu1_covered95),
    Pi = mean(Pi_covered95)
    ) %>%
  arrange(Prior, size_tag, miss_tag)

fact <- summ$Prior
summ <- summ %>% ungroup %>%
  select(-Prior) %>%
  rename(
    Missing = miss_tag,
    Size    = size_tag
    ) %>%
  split(., fact) %>%
  `attr<-`("subheadings", paste(names(.), "Prior")) %>%
  xtableList(
    caption = "Coverage probability at the 95\\% level by prior (right/wrong), size of the tree, and proportion of missingness.  Estimations with the \\emph{right} prior use the same priors as the data generating process, whereas estimations with the \\emph{wrong} prior used a prior that had a mean twice as large as the data generating process.",
    label   = "tab:coverage95-method-size-missigness") %>%
  print %>%
  cat(file = "tables/03_coverage95_by_method_size_missingness.tex")


summ <- bias %>% 
  group_by(Prior, miss_tag) %>%
  summarize(
    psi0 = mean(psi0_covered95),
    psi1 = mean(psi1_covered95),
    mu0 = mean(mu0_covered95),
    m1 = mean(mu1_covered95),
    Pi = mean(Pi_covered95)
  ) %>%
  arrange(Prior, miss_tag)

fact <- summ$Prior
summ %>% ungroup %>%
  select(-Prior) %>%
  rename(Missing = miss_tag) %>%
  split(., fact) %>%
  `attr<-`("subheadings", paste(names(.), "Prior")) %>%
  xtableList(
    caption = "Coverage probability at the 95\\% level by prior (right/wrong), and proportion of missingness. Estimations with the \\emph{right} prior use the same priors as the data generating process, whereas estimations with the \\emph{wrong} prior used a prior that had a mean twice as large as the data generating process.",
    label   = "tab:coverage95-method-missigness") %>%
  print %>%
  cat(file = "tables/03_coverage95_by_method_missingness.tex")

