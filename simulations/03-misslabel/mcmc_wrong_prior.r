#!/bin/sh
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=03-mcmc_wrong_prior

library(aphylo)
library(slurmR)
source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds") #[1:NSAMPLES]
dat0 <- lapply(dat0, "[[", "atree_m")

# Setting the seed
set.seed(111222)

parnames <- c("psi0", "psi1", "mu_d0", "mu_d1", "mu_s0", "mu_s1", "Pi")
mcmc.par   <- matrix(
  runif(length(parnames) * mcmc.nchains),
  ncol=length(parnames)
)
mcmc.prior <- function(p) {
  dbeta(p, ALPHA_PAR_WRONG[parnames], BETA_PAR[parnames])
}

job <- Slurm_lapply(
  dat0,
  mcmc_lite,
  model      = dat ~ psi + mu_d + mu_s + Pi,
  params     = mcmc.par,
  priors     = mcmc.prior,
  nsteps     = mcmc.nsteps,
  nchains    = mcmc.nchains,
  burnin     = mcmc.burnin,
  thin       = mcmc.thin,
  njobs      = NJOBS,
  mc.cores   = 1L,
  multicore  = mcmc.multicore, # TRUE,
  job_name   = "03-mcmc_wrong_prior-lapply",
  plan       = "submit",
  export     = ls()
)

saveRDS(job, paste0(PROJECT_PATH, "/simulations/03-misslabel/mcmc_wrong_prior-job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/03-misslabel/mcmc_wrong_prior.rds"
    )
  )

job
