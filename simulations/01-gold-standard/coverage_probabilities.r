rm(list = ls())

library(dplyr)
library(tidyselect)
library(magrittr)
library(xtable)

# Setting xtable options
options(xtable.booktabs = TRUE)
options(xtable.include.rownames = FALSE)
options(xtable.sanitize.text.function = function(x) x)

bias <- readRDS("simulations/01-gold-standard/bias.rds")

# Coverage probability by size -------------------------------------------------
summ <- bias %>% 
  group_by(Prior, size_tag) %>%
  summarize_at(vars(ends_with("covered95")), ~ mean(.)) %>%
  set_colnames(gsub("[_]covered95$", "", colnames(.))) %>%
  arrange(desc(Prior), size_tag) %>%
  # Cleaning latex errors
  # See https://www.sharelatex.com/learn/Errors/Illegal_unit_of_measure_(pt_inserted)
  ungroup %>%
  mutate(Size = paste0("{", size_tag, "}")) %>%
  select(-size_tag) %>%
  rename(
    `$\\mu_{d01}$` = mu_d0,
    `$\\mu_{d10}$` = mu_d1,
    `$\\mu_{s01}$` = mu_s0,
    `$\\mu_{s10}$` = mu_s1,
    # `$\\psi_{01}$` = psi0,
    # `$\\psi_{10}$` = psi1,
    `$\\pi$` = Pi,
  )

fact <- summ$Prior
summ %>% 
  select(-Prior) %>%
  split(., fact) %>%
  `attr<-`("subheadings", paste(names(.), "Prior")) %>%
  xtableList(
    caption = paste(
      "Coverage probability at the 95\\% level by prior (right/wrong) and size",
      "of the tree for the fully annotated model."
    ),
    label   = "tab:1coverage95-prior-size") %>%
  print %>%
  cat(file = "tables/01_coverage95_by_prior_and_size.tex")

