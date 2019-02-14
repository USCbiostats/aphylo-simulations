library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds") #[1:NSAMPLES]
dat0 <- lapply(dat0, "[[", "atree_m")

# Setting the seed
set.seed(111222)

mcmc.par   <- matrix(runif(5*mcmc.nchains), ncol=5)
mcmc.prior <- function(p) {
  dbeta(p, c(4, 4, 4, 4, 4), c(16, 16, 16, 16, 16))
}

job <- Slurm_lapply(
  dat0,
  mcmc_lite,
  model      = dat ~ psi + mu + Pi,
  params     = mcmc.par,
  priors     = mcmc.prior,
  nbatch     = mcmc.nbatch,
  nchains    = mcmc.nchains,
  burnin     = mcmc.burnin,
  thin       = mcmc.thin,
  njobs      = 55L,
  mc.cores   = 4L,
  multicore  = mcmc.multicore, # TRUE,
  job_name   = "02-missing-wrong-prior",
  submit     = TRUE
)

saveRDS(job, paste0(PROJECT_PATH, "/simulations/02-missing/wrong-prior-estimates-job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/02-missing/mcmc_wrong_prior_estimates.rds"
    )
  )

res
job
