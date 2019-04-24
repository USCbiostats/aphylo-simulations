# Local paths
source("global-paths.r")
# STAGING_PATH <- "/staging/pdt/vegayon/aphylo-simulations"
# PANTHER_PATH <- "/auto/pmd-02/pdt/pdthomas/panther/famlib/rel/PANTHER13.1_altVersion/hmmscoring/PANTHER13.1/books"
# PROJECT_PATH <- "/home/rcf-proj2/pdt/vegayon/aphylo-simulations"

opts_sluRm$set_chdir(STAGING_PATH)
opts_sluRm$set_opts(
  partition     = "thomas",
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
ALPHA_PAR <- c( 2, 2, 2, 2, 7,38, 2)
BETA_PAR  <- c(38,38,38,38, 3, 2,38)

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
      autostop     = 5e3
    )
    ,
    priors            = priors,
    check.informative = FALSE,
    reduced_pseq      = reduced_pseq.
  ),
  error = function(e) e
  )
  
  # If it didn't worked, then return error
  if (inherits(ans, "error")) return(ans)
  
  # ans[["dat"]] <- NULL
  
  return(ans)
  
  
}

