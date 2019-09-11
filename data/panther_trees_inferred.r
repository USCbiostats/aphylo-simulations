# library(dplyr)
# library(tidyr)
# library(magrittr)
library(aphylo)
library(parallel)

# source("global-paths.r")

# Reading the list of candidate functions and keeping the trees only
candidate_functions <- data.table::fread("data/candidate_functions_inferred.csv")
candidate_trees     <- unique(candidate_functions$substring)

# Reading panther trees
trees <- mclapply(
  sprintf("%s/%s/tree.tree", PANTHER_PATH, candidate_trees), 
  read_panther,
  mc.cores   = 4 #,
  # job_name   = "candidate-trees",
  # tmp_path   = STAGING_PATH,
  # sbatch_opt = list(account = "lc_pdt", partition = "thomas"),
  # plan       = "wait"
)

# trees <- Slurm_collect(trees)

saveRDS(trees, "data/panther_trees_inferred.rds", compress = TRUE)
