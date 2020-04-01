#!/bin/sh/
#SBATCH --account=lc_pdt
#SBATCH --partition=thomas
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=one-tree-at-a-time
#SBATCH --output=one-tree-at-a-time.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=g.vegayon@gmail.com

library(slurmR)
library(aphylo)

cl <- makeSlurmCluster(
  40L,
  account   = "lc_pdt",
  partition = "thomas",
  time      = "04:00:00",
  `mem-per-cpu` = "2G"
)

# Loading the aphylo package
setup <- tryCatch(clusterEvalQ(cl, {
  library(slurmR)
  library(aphylo)
  library(coda)
  
  shrink_towards_half <- function(x, margin=.01) {
    
    x[x < (.5 - margin)] <- x[x < (.5 - margin)] + margin
    x[x > (.5 + margin)] <- x[x > (.5 + margin)] - margin
    
    x
  }
  
  # Common parameters
  prior.  <- bprior(c(2,2,9,9,2,2,2,2,2), c(9,9,2,2,9,9,9,9,9))
  warmup. <- 1000
  freq.   <- 10
  lb.     <- 1e-5
  ub.     <- 1 - 1e-5
  
  mcmc. <- list(
    nchains      = 4L,
    multicore    = FALSE, 
    burnin       =  5000L,
    nsteps       = 10000L,
    conv_checker = NULL,
    kernel       = fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup.),
    thin         = 10L
  )
}), error = function(e) e)

if (inherits(setup, "error")) {
  stopCluster(cl)
  stop("An error ocurred during setup. ", setup)
}


set.seed(1362)

# Case 1: (Almost) Fully annotated trees ---------------------------------------

# Data preprocessing
trees <- readRDS("../data/candidate_trees.rds")
trees <- do.call(c, trees)
trees <- unlist(lapply(trees, function(d) {
  lapply(1:Nann(d), function(i) d[,i])
}), recursive = FALSE)
trees <- do.call(c, trees)

# # Selecting using balance
# b <- balance_ann(trees)
# trees <- trees[(b >= quantile(b, .2))]

# 1.A: No prior
res <- parLapply(
  cl, trees, function(t.) {
    
    # Finding MLEs
    tryCatch({
      mle  <- aphylo_mle(t. ~ psi + mu_d + mu_s + Pi, priors = prior.)
      
      # Estimating MCMC
      mcmc.$kernel <- fmcmc::kernel_adapt(
        lb = lb., ub = ub., warmup = warmup. #, freq = freq.
        )
      
      mcmc <- aphylo_mcmc(
        t. ~ psi + mu_d + mu_s + Pi,
        priors  = prior.,
        control = mcmc.,
        params  = shrink_towards_half(coef(mle))
      )
      
      # Returning prediction
      list(
        ps = prediction_score(mcmc, loo = TRUE),
        coef = coef(mcmc),
        vcov = vcov(mcmc),
        gelman = coda::gelman.diag(mcmc$hist)
      )
    }, error = function(e) e)
    
  }
)


saveRDS(res, "estimates_disjoint.rds")
stopCluster(cl)

