# library(dplyr)
# library(tidyr)
# library(magrittr)
library(aphylo)
library(slurmR)

source("global-paths.r")

# Reading the list of candidate functions and keeping the trees only
candidate_functions <- readr::read_csv("data/candidate_functions.csv")
candidate_trees     <- unique(candidate_functions$substring)

# Reading panther trees
trees <- Slurm_lapply(
  sprintf("%s/%s/tree.tree", PANTHER_PATH, candidate_trees), 
  read_panther,
  njobs      = 40,
  job_name   = "candidate-trees",
  tmp_path   = STAGING_PATH,
  sbatch_opt = list(account = "lc_pdt", partition = "thomas"),
  plan       = "wait"
)

trees <- Slurm_collect(trees)

saveRDS(trees, "data/panther_trees.rds", compress = TRUE)
