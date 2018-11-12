library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds") #:NSAMPLES]
dat0 <- lapply(dat0, "[[", "atree_o")

# Setting the seed
set.seed(111222)

mcmc.par   <- matrix(runif(7*mcmc.nchains), ncol=7)
mcmc.prior <- function(p) {
  dbeta(p, c(2, 2, 2, 2, 8, 38, 2), c(18, 18, 18, 18, 2, 4, 18))
}


job <- Slurm_lapply(
  dat0,
  mcmc_lite,
  model      = dat ~ psi + eta + mu + Pi,
  params     = mcmc.par,
  priors     = mcmc.prior,
  nsteps     = mcmc.nsteps,
  nchains    = mcmc.nchains,
  burnin     = mcmc.burnin,
  thin       = mcmc.thin,
  reduced_pseq. = FALSE,
  njobs      = 55L,
  mc.cores   = 4L,
  multicore  = mcmc.multicore, # TRUE,
  job_name   = "04-full-model-wrong-prior",
  job_path   = STAGING_PATH,
  submit     = TRUE
)

saveRDS(job, paste0(PROJECT_PATH, "/simulations/04-full-model/wrong-prior-estimates-job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/04-full-model/mcmc_wrong_prior_estimates.rds"
    )
  )

# res
# job
