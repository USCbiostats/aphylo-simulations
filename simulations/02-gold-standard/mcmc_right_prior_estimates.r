library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat <- readRDS("simulations/dgp.rds")

# Priors and starting point
mcmc.par   <- c(0.1, 0.1, 0.1, 0.1, 0.7, 0.9, 0.1)
mcmc.prior <- function(p) {
  dbeta(p, c(2, 2, 2, 2, 7, 18, 2), c(18, 18, 18, 18, 3, 2, 18))
}

# Setting the seed
set.seed(1)

job <- Slurm_lapply(
    dat,
    mcmc_lite,
    par      = mcmc.par,
    priors   = mcmc.prior,
    nbatch   = mcmc.nbatch,
    nchains  = mcmc.nchains,
    burnin   = mcmc.burnin,
    thin     = mcmc.thin,
    njobs    = 9,
    mc.cores = 10L,
    multicore= mcmc.multicore, # TRUE,
    job_name = "mcmc_right_prior",
    job_path = STAGING_PATH,
    submit = TRUE,
    sbatch_opt = list(ntasks = 10, time="02:00:00", `cpus-per-task` = 1)
  )

saveRDS(job, paste0(PROJECT_PATH, "/simulations/02-gold-standard/job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/02-gold-standard/mcmc_right_prior_estimates.rds"
    )
  )

res
job
