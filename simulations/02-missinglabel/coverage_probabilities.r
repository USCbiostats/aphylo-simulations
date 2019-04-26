rm(list = ls())

library(dplyr)
library(tidyselect)
library(magrittr)
library(xtable)

# Setting xtable options
options(xtable.booktabs = TRUE)
options(xtable.include.rownames = FALSE)
options(xtable.sanitize.text.function = function(x) x)

bias <- readRDS("simulations/02-missinglabel/bias.rds")

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
    `$\\mu_{01}$` = mu0,
    `$\\mu_{10}$` = mu1,
    `$\\psi_{01}$` = psi0,
    `$\\psi_{10}$` = psi1,
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
      "of the tree for the partially annotated model."
    ),
    label   = "tab:2coverage95-prior-size") %>%
  print %>%
  cat(file = "tables/02_coverage95_by_prior_and_size.tex")


# Coverage probability by size -------------------------------------------------
summ <- bias %>% 
  group_by(Prior, miss_tag) %>%
  summarize_at(vars(ends_with("covered95")), ~ mean(.)) %>%
  set_colnames(gsub("[_]covered95$", "", colnames(.))) %>%
  arrange(desc(Prior), miss_tag) %>%
  # Cleaning latex errors
  # See https://www.sharelatex.com/learn/Errors/Illegal_unit_of_measure_(pt_inserted)
  ungroup %>%
  mutate(Missing = paste0("{", miss_tag, "}")) %>%
  select(-miss_tag) %>%
  rename(
    `$\\mu_{01}$` = mu0,
    `$\\mu_{10}$` = mu1,
    `$\\psi_{01}$` = psi0,
    `$\\psi_{10}$` = psi1,
    `$\\pi$` = Pi,
  )

fact <- summ$Prior
summ %>% 
  select(-Prior) %>%
  split(., fact) %>%
  `attr<-`("subheadings", paste(names(.), "Prior")) %>%
  xtableList(
    caption = paste(
      "Coverage probability at the 95\\% level by prior (right/wrong) and",
      "proportion of missing labels for the partially annotated model."
    ),
    label   = "tab:2coverage95-prior-missingness") %>%
  print %>%
  cat(file = "tables/02_coverage95_by_prior_and_missingness.tex")

