salloc --mem=4G --time=02:00:00 --mail-type=ALL \
	R CMD BATCH --vanilla simulations/02-missing/mcmc_right_prior_estimates.r &
