salloc --mem=4G --time=02:00:00 --mail-type=ALL --quiet \
	R CMD BATCH --vanilla simulations/04-full-model/mcmc_right_prior_estimates.r \
	simulations/04-full-model/mcmc_right_prior_estimates.rout &
