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

# Filtering for candidate functions
candidate_functions <- counts %>%
  filter(no > 0, yes > 0) %>%
  arrange(desc(no))

# Saving
readr::write_csv(candidate_functions, path = "data/candidate_functions.csv")
