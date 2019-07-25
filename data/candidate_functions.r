library(dplyr)
library(magrittr)

# Reading true data
true_annotations <- readr::read_delim(
  file = "data-raw/true_annotations",
  delim=";", col_names = TRUE
  )

# Counting how many functional annotations are present per tree, per function.
counts <- true_annotations %>%
  mutate(present = if_else(
    is.na(qualifier), TRUE,
    if_else(qualifier == "NOT", FALSE, TRUE))
    ) %>%
  group_by(substring, term) %>%
  summarize(
    yes = sum(present),
    no  = n() - yes
  )

# General characteristics of the panther dataset -------------------------------
candidate_functions_stats <- counts %>%
  ungroup %>%
  group_by(substring) %>%
  mutate(n_annotations = n()) %>%
  filter(1:n() == 1) %>%
  ungroup() %$%
  table(n_annotations) %>%
  as_tibble %>%
  mutate(
    n_annotations = as.integer(n_annotations),
    flag          = n_annotations >= 10
  ) %>%
  arrange(n_annotations) %>%
  group_by(flag) %>%
  mutate(tot = sum(n)*flag) %>%
  ungroup %>%
  mutate(
    n      = if_else(1:n() == 10, tot, n),
    naccum = cumsum(n)
  ) %>%
  select(n_annotations, n, naccum) %>%
  filter(1:n() <= 10) %>%
  mutate(
    naccum_pcent  = naccum / naccum[n()],
    n_annotations = as.character(n_annotations),
    n_annotations = if_else(n_annotations == "10", "10 or more", n_annotations)
  ) %>%
  rename(
    `N Ann. per Tree` = n_annotations,
    `N`               = n,
    `Cum. count`      = naccum,
    `Cum. prop.`      = naccum_pcent
  )

xtable::xtable(candidate_functions_stats) %>%
  xtable::`caption<-`(paste(
    "Distribution of number of functions per tree in PantherDB."
  )) %>%
  xtable::print.xtable(
    file = "data/candidate_functions_stats.tex",
    include.rownames = FALSE,
    booktabs = TRUE
  )

# Filtering for candidate functions
candidate_functions <- counts %>%
  filter(no > 0, yes > 0) %>%
  arrange(desc(no))

# Saving
readr::write_csv(candidate_functions, path = "data/candidate_functions.csv")
