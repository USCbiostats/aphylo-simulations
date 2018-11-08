library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds")[1:NSAMPLES]
dat0 <- lapply(dat0, "[[", "atree")

# Setting the seed
set.seed(111222)

mcmc.par   <- matrix(runif(7*mcmc.nchains), ncol=7)
mcmc.prior <- function(p) {
  dbeta(p, c(2, 2, 2, 2, 7, 18, 2), c(18, 18, 18, 18, 3, 2, 18))
}


job <- Slurm_lapply(
  dat0,
  mcmc_lite,
  model      = dat ~ psi + eta + mu + Pi,
  params     = mcmc.par,
  priors     = mcmc.prior,
  nbatch     = mcmc.nbatch,
  nchains    = mcmc.nchains,
  burnin     = mcmc.burnin,
  thin       = mcmc.thin,
  njobs      = 40,
  mc.cores   = 5L,
  multicore  = mcmc.multicore, # TRUE,
  job_name   = "01-gold-standard-right-prior",
  submit     = TRUE
)

saveRDS(job, paste0(PROJECT_PATH, "/simulations/01-gold-standard/job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/01-gold-standard/mcmc_right_prior_estimates.rds"
    )
  )

head(res)

job
