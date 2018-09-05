library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat <- readRDS("simulations/dgp.rds")[1:100]

# Priors and starting point
mcmc.par   <- c(0.1, 0.1, 0.1, 0.1, 0.7, 0.9, 0.1)
mcmc.prior <- function(p) {
  dbeta(p, c(2, 2, 2, 2, 7, 18, 2), c(18, 18, 18, 18, 3, 2, 18))
}

job <- Slurm_lapply(
    dat,
    mcmc_lite,
    par      = mcmc.par,
    priors   = mcmc.prior,
    nbatch   = mcmc.nbatch,
    nchains  = mcmc.nchains,
    burnin   = mcmc.burnin,
    thin     = mcmc.thin,
    njobs    = 10,
    mc.cores = 8L,
    multicore= mcmc.multicore,
    job_name = "mcmc_right_prior2",
    job_path = paste0(PROJECT_PATH, "/simulations/02-gold-standard/"),
    submit = TRUE
  )

saveRDS(Slurm_collect(job), "simulations/02-gold-standard/mcmc_right_prior_estimates.rds")

