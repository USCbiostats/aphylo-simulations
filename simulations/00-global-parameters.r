# Local paths
source("global-paths.r")

opts_slurmR$set_tmp_path(STAGING_PATH)
opts_slurmR$set_opts(
  partition     = "scavenge",
  account       = "lc_pdt",
  time          = "05:00:00",
  `mem-per-cpu` = "2G"
)

if (!dir.exists(STAGING_PATH))
  dir.create(STAGING_PATH)

NSAMPLES     <- 10000

# MCMC
mcmc.nsteps  <- 1e5
mcmc.burnin  <- 2e4
mcmc.thin    <- 50
mcmc.nchains <- 4
mcmc.multicore <- FALSE

# True DGP parameters
# These are given in the following way
# psi0, psi1, mu_d0, mu_d1, mu_s0, mu_s1, eta0, eta1, Pi
ALPHA_PAR <- c( 2,  2, 38, 10,  2,  2, 38, 38,  2)
BETA_PAR  <- c(38, 38,  2, 10, 38, 38,  2,  2, 38)

# Function to estimate model using MCMC
mcmc_lite <- function(
  dat,
  model,
  params,
  priors  = NULL,
  nsteps  = mcmc.nsteps,
  nchains = mcmc.nchains,
  burnin  = mcmc.burnin,
  thin    = mcmc.thin,
  multicore = mcmc.multicore,
  reduced_pseq. = TRUE
) {
  
  # Making sure it runs in the same place
  environment(model) <- environment()

   # Try to estimate the model
  ans <- tryCatch(aphylo_mcmc(
    model,
    params  = params,
    control = list(
      nsteps       = nsteps,
      nchains      = nchains,
      burnin       = burnin,
      thin         = thin,
      multicore    = multicore,
      conv_checker = fmcmc::convergence_gelman(1.05),
      autostop     = 1e4
    )
    ,
    priors            = priors,
    check_informative = FALSE,
    reduced_pseq      = reduced_pseq.
  ),
  error = function(e) e
  )
  
  # If it didn't worked, then return error
  if (inherits(ans, "error")) return(ans)
  
  # ans[["dat"]] <- NULL
  
  return(ans)
  
  
}

