library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds")
dat0 <- lapply(dat0, "[[", "atree_o")

# Priors and starting point
mcmc.par   <- c(0.1, 0.1, 0.1, 0.1, 0.1)
mcmc.prior <- function(p) {
  dbeta(p, c(2, 2, 2, 2, 2), c(18, 18, 18, 18, 18))
}

# Setting the seed
set.seed(1)

job <- Slurm_lapply(
    dat0,
    mcmc_lite,
    model      = dat ~ mu + psi + Pi,
    par        = mcmc.par,
    priors     = mcmc.prior,
    nbatch     = mcmc.nbatch,
    nchains    = mcmc.nchains,
    burnin     = mcmc.burnin,
    thin       = mcmc.thin,
    njobs      = 9,
    mc.cores   = 10L,
    multicore  = mcmc.multicore, # TRUE,
    job_name   = "aphylo-03-pub-bias",
    job_path   = STAGING_PATH,
    submit     = TRUE,
    sbatch_opt = list(
      ntasks          =  10,
      time            = "02:00:00",
      `cpus-per-task` = 1
      )
  )

saveRDS(job, paste0(PROJECT_PATH, "/simulations/03-pub-bias/right-prior-estimates-job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/03-pub-bias/mcmc_right_prior_estimates.rds"
    )
  )

# res
# job
