library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds") #[1:NSAMPLES]
dat0 <- lapply(dat0, "[[", "atree_i")

# Setting the seed
set.seed(111222)

mcmc.par   <- matrix(runif(5*mcmc.nchains), ncol=5)
mcmc.prior <- function(p) {
  dbeta(p, c(2, 2, 2, 2, 2), c(38, 38, 38, 38, 38))
}

job <- Slurm_lapply(
  dat0,
  mcmc_lite,
  model      = dat ~ psi + mu + Pi,
  params     = mcmc.par,
  priors     = mcmc.prior,
  nsteps     = mcmc.nsteps,
  nchains    = mcmc.nchains,
  burnin     = mcmc.burnin,
  thin       = mcmc.thin,
  njobs      = 55L,
  mc.cores   = 4L,
  multicore  = mcmc.multicore, # TRUE,
  job_name   = "02-missinglabel-right-prior",
  submit     = TRUE
)

saveRDS(job, paste0(PROJECT_PATH, "/simulations/02-missinglabel/right-prior-estimates-job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/02-missinglabel/mcmc_right_prior_estimates.rds"
    )
  )

res
job
