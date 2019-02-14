library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds") #:NSAMPLES]
dat0 <- lapply(dat0, "[[", "atree_o")

# Setting the seed
set.seed(111222)

mcmc.par   <- matrix(runif(7*mcmc.nchains), ncol=7)
mcmc.prior <- function(p) {
  dbeta(p, c(2, 2, 2, 2, 7, 38, 2), c(38, 38, 38, 38, 3, 2, 38))
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
  job_name   = "04-full-model-right-prior",
  job_path   = STAGING_PATH,
  submit     = TRUE
)

saveRDS(job, paste0(PROJECT_PATH, "/simulations/04-full-model/right-prior-estimates-job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/04-full-model/mcmc_right_prior_estimates.rds"
    )
  )

# res
# job
