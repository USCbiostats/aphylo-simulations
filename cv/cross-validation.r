library(aphylo)
library(sluRm)

# source("global-paths.r")
source("simulations/00-global-parameters.r")
trees <- readRDS("data/candidate_trees.rds")

# Unlisting functions
trees <- unlist(lapply(trees, function(d) lapply(1:Nann(d), function(i) d[i])), recursive = FALSE)

f <- function(d) {
  
  model <- d ~ psi + Pi + mu
  environment(model) <- environment()
  
  tryCatch(aphylo_cv(
    model,
    control = list(
      nsteps       = 1e5L,
      nchains      = 4L,
      burnin       = 2e4L,
      thin         = 50L,
      multicore    = FALSE,
      conv_checker = amcmc::gelman_convergence(1.05),
      autostop     = 5e3L
    ),
    priors = bprior(2,20),
    reduced_pseq = TRUE
    ), error = function(e) e)
  
}


ans_slurm <- Slurm_lapply(trees, f, njobs=100, mc.cores=4L, job_name="aphylo-cross-validation")
saveRDS(ans_slurm, "cv/cross_validation-job.rds")
saveRDS(Slurm_collect(ans_slurm), "cv/cross_validation.rds")



