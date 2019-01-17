
source("global-paths.r")
candidate_functions <- readr::read_csv("data/candidate_functions.csv")
candidate_trees     <- unique(candidate_functions$substring)

library(aphylo)

trees <- lapply(
  sprintf("%s/%s/tree.tree", PANTHER_PATH, candidate_trees[1:10]),
  read_panther
  )
