library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds")[1:NSAMPLES]
dat0 <- lapply(dat0, "[[", "atree")

# Setting the seed
set.seed(111222)

mcmc.par   <- matrix(runif(7*mcmc.nchains), ncol=7)
mcmc.prior <- function(p) {
  dbeta(p, c(4, 4, 4, 4, 7, 18, 4), c(18, 18, 18, 18, 1, 1, 18))
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
  njobs      = 15,
  mc.cores   = 10L,
  multicore  = mcmc.multicore, # TRUE,
  job_name   = "01-gold-standard-wrong-prior",
  job_path   = STAGING_PATH,
  submit     = TRUE,
  sbatch_opt = list(
    ntasks          =  10L,
    time            = "02:00:00",
    `cpus-per-task` = 1,
    `mem-per-cpu`   = "500M",
    account         = "lc_pdt",
    partition       = "thomas"
    )
)

saveRDS(job, paste0(PROJECT_PATH, "/simulations/01-gold-standard/job-wrong-prior.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/01-gold-standard/mcmc_wrong_prior_estimates.rds"
    )
  )

head(res)

job
