rm(list = ls())
library(aphylo)
load("simulations/data_and_functions.rda")

if (file.exists("simulations/mcmc_right_prior_estimates.rda")) {
  load("simulations/mcmc_right_prior_estimates.rda")
  start <- i - 1L
  set.seed(curseed_MCMC_right_prior)
} else {
  set.seed(1223)
  ans_MCMC_right_prior <- vector("list", nsim)
  start <- 1
}

# Priors and starting point
mcmc.par   <- c(rep(1/40, 2), rep(1/20, 3))
mcmc.prior <- function(p) {
  c(dbeta(p[1:2], 1, 39), dbeta(p[3:5], 1, 19))
}

for (i in start:nsim) {
  
  # MCMC estimators
  ans_MCMC_right_prior[[i]] <- mcmc_lite(
    dat    = dat_obs[[i]], 
    par    = mcmc.par,
    priors = mcmc.prior
    )
  
  # Printing on screen (and saving)
  if (!(i %% 10)) {
    message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))
    
    # Just in case we need to restart this... so we can continue where we were
    curseed_MCMC_right_prior <- .Random.seed
    
    # Storing all objects
    save(
      curseed_MCMC_right_prior, ans_MCMC_right_prior, i,
      file     = "simulations/mcmc_right_prior_estimates.rda",
      compress = FALSE
      )
    
  } else message(".", appendLF = FALSE)
   
}

# Terminating
message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))

curseed_MCMC_right_prior <- .Random.seed
save(
  curseed_MCMC_right_prior, ans_MCMC_right_prior, i,
  file     = "simulations/mcmc_right_prior_estimates.rda",
  compress = FALSE
  )
