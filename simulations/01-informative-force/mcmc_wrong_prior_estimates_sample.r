set.seed(18888)

load("simulations/mcmc_wrong_prior_estimates.rda")

N   <- 1e3
mcmc_wrong_prior_estimates_sample <- ans_MCMC_wrong_prior[1:N]

save(mcmc_wrong_prior_estimates_sample,
    file = "simulations/mcmc_wrong_prior_estimates_sample.rda")
