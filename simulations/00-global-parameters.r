# Local paths
PANTHER_PATH <- "/auto/pmd-02/pdt/pdthomas/panther/famlib/rel/PANTHER13.1_altVersion/hmmscoring/PANTHER13.1/books"
PROJECT_PATH <- "/home/rcf-proj2/pdt/vegayon/aphylo-simulations"

NSAMPLES     <- 1e4

# MCMC
mcmc.nbatch  <- 2e4
mcmc.burnin  <- 5e3
mcmc.thin    <- 20
mcmc.nchains <- 4
mcmc.parallel <- FALSE

# True DGP parameters
ALPHA_PAR <- c( 2, 2, 2, 2, 7,18, 2)
BETA_PAR  <- c(18,18,18,18, 3, 2,18)

# Function to estimate model using MCMC
mcmc_lite <- function(
  dat,
  par,
  priors  = NULL,
  nbatch  = mcmc.nbatch,
  nchains = mcmc.nchains,
  burnin  = mcmc.burnin,
  thin    = mcmc.thin,
  parallel = mcmc.parallel
) {
  
  cl <- parallel::makeForkCluster(nchains)
  
  # Try to estimate the model
  ans <- tryCatch(aphylo_mcmc(
    dat$atree ~ psi + mu + eta + Pi,
    control = list(
      nbatch = nbatch, nchains=nchains, cl = cl, burnin = burnin, thin=thin,
      parallel = mcmc.parallel),
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

