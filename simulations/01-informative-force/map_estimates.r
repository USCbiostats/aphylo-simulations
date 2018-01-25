rm(list = ls())

load("simulations/data_and_functions.rda")
library(aphylo)

set.seed(1223)
ans_MAP <- vector("list", nsim)

# Priors and starting point
mcmc.par   <- c(rep(1/40, 2), rep(1/20, 3))
mcmc.prior <- function(p) {
  c(dbeta(p[1:2], 1, 39), dbeta(p[3:5], 1, 19))
}


for (i in 1:nsim) {
  
  # MAP estimators
  ans_MAP[[i]] <- mle_lite(dat_obs[[i]], FALSE, priors = mcmc.prior, par = mcmc.par)
  
  # Printing on screen (and saving)
  if (!(i %% 10)) {
    message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))
    
    # Just in case we need to restart this... so we can continue where we were
    curseed_MAP <- .Random.seed
    
    # Storing all objects
    save(curseed_MAP, ans_MAP, i, file = "simulations/map_estimates.rda", compress = FALSE)
  } else message(".", appendLF = FALSE)
  
}

# Terminating
message(sprintf("Simulation %04i/%04i (%6.2f %%) complete.", i, nsim, i/nsim*100))

curseed_MAP <- .Random.seed
save(curseed_MAP, ans_MAP, i, file = "simulations/map_estimates.rda", compress = FALSE)
