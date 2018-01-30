rm(list = ls())
library(aphylo)
load("simulations/03_uniform_drop_of_0/data_and_functions.rda")

if (file.exists("simulations/03_uniform_drop_of_0/mcmc_wrong_prior_estimates.rda")) {
  load("simulations/03_uniform_drop_of_0/mcmc_wrong_prior_estimates.rda")
  start <- i - 1L
  set.seed(curseed_MCMC_wrong_prior)
} else {
  set.seed(1223)
  ans_MCMC_wrong_prior <- vector("list", nsim)
  start <- 1
}

# Priors and starting point
mcmc.par   <- c(rep(2/40, 2), rep(2/20, 3))
mcmc.prior <- function(p) {
  c(dbeta(p[1:2], 2, 39), dbeta(p[3:5], 2, 19))
}

for (i in 1:nsim) {
  
  # MCMC estimators
  ans_MCMC_wrong_prior[[i]] <- mcmc_lite(
    dat_obs[[i]], mcmc.par, priors = mcmc.prior
  )
  
  # Printing on screen (and saving)
  if (!(i %% 10)) {
    message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))
    
    # Just in case we need to restart this... so we can continue where we were
    curseed_MCMC_wrong_prior <- .Random.seed
    
    # Storing all objects
    save(
      curseed_MCMC_wrong_prior, ans_MCMC_wrong_prior, i,
      file = "simulations/03_uniform_drop_of_0/mcmc_wrong_prior_estimates.rda",
      compress = FALSE
      )
  } else message(".", appendLF = FALSE)
  
}

# Terminating
message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))

curseed_MCMC_wrong_prior <- .Random.seed
save(curseed_MCMC_wrong_prior, ans_MCMC_wrong_prior, i,
     file = "simulations/03_uniform_drop_of_0/mcmc_wrong_prior_estimates.rda",
     compress = FALSE
     )
