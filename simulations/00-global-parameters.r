# Local paths
STAGING_PATH <- "/staging/pdt/vegayon/aphylo-simulations"
PANTHER_PATH <- "/auto/pmd-02/pdt/pdthomas/panther/famlib/rel/PANTHER13.1_altVersion/hmmscoring/PANTHER13.1/books"
PROJECT_PATH <- "/home/rcf-proj2/pdt/vegayon/aphylo-simulations"

if (!dir.exists(STAGING_PATH))
  dir.create(STAGING_PATH)

NSAMPLES     <- 100

# MCMC
mcmc.nbatch  <- 5e4
mcmc.burnin  <- 1e4
mcmc.thin    <- 100
mcmc.nchains <- 2
mcmc.multicore <- FALSE

# True DGP parameters
ALPHA_PAR <- c( 2, 2, 2, 2, 7,18, 2)
BETA_PAR  <- c(18,18,18,18, 3, 2,18)

# Function to estimate model using MCMC
mcmc_lite <- function(
  dat,
  model,
  par,
  priors  = NULL,
  nbatch  = mcmc.nbatch,
  nchains = mcmc.nchains,
  burnin  = mcmc.burnin,
  thin    = mcmc.thin,
  multicore = mcmc.multicore
) {
  
   # Try to estimate the model
  ans <- tryCatch(aphylo_mcmc(
    model,
    control = list(
      nbatch    = nbatch,
      nchains   = nchains,
      burnin    = burnin,
      thin      = thin,
      multicore = multicore
    )
    ,
    priors = priors,
    check.informative = FALSE
  ),
  error = function(e) e
  )
  
  # If it didn't worked, then return error
  if (inherits(ans, "error")) return(ans)
  
  # ans[["dat"]] <- NULL
  
  return(ans)
  
  
}

