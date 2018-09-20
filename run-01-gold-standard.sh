salloc --mem=24G --time=12:00:00 --mail-type=ALL --quiet --job-name=01-gold-standard \
	R CMD BATCH --vanilla simulations/01-gold-standard/mcmc_right_prior_estimates.r \
	simulations/01-gold-standard/mcmc_right_prior_estimates.rout &
