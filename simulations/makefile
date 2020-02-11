slurmr := ~/./slurmr
.PHONY: all
01:=01-gold-standard
02:=02-missinglabel
03:=03-misslabel
globalpars:=00-global-parameters.r

all: */mcmc_right_prior.rds */mcmc_wrong_prior.rds 
01right: $(01)/mcmc_right_prior.rds

# Data generating process
dgp.rds: 01-dgp.r $(globalpars)
	$(slurmr) 01-dgp.r

# Running the gold standard ----------------------------------------------------
$(01)/mcmc_right_prior.rds: $(01)/mcmc_right_prior.r dgp.rds $(globalpars)
	$(slurmr) $(01)/mcmc_right_prior.r
$(01)/mcmc_wrong_prior.rds: $(01)/mcmc_wrong_prior.r dgp.rds $(globalpars)
	$(slurmr) $(01)/mcmc_wrong_prior.r
$(01)/bias.rds: $(01)/bias.r $(01)/mcmc_wrong_prior.rds $(01)/mcmc_right_prior.rds 00-bias-functions.r
	$(slurmr) $(01)/bias.r

# Running simulations with missing data ----------------------------------------
$(02)/mcmc_right_prior.rds: $(02)/mcmc_right_prior.r dgp.rds $(globalpars)
	$(slurmr) $(02)/mcmc_right_prior.r
$(02)/mcmc_wrong_prior.rds: $(02)/mcmc_wrong_prior.r dgp.rds $(globalpars)
	$(slurmr) $(02)/mcmc_wrong_prior.r
$(02)/bias.rds: $(02)/bias.r $(02)/mcmc_wrong_prior.rds $(02)/mcmc_right_prior.rds 00-bias-functions.r
	$(slurmr) $(02)/bias.r


# Running simulations with missing data and mis labels -------------------------
$(03)/mcmc_right_prior.rds: $(03)/mcmc_right_prior.r dgp.rds $(globalpars)
	$(slurmr) $(03)/mcmc_right_prior.r
$(03)/mcmc_wrong_prior.rds: $(03)/mcmc_wrong_prior.r dgp.rds $(globalpars)
	$(slurmr) $(03)/mcmc_wrong_prior.r
$(03)/bias.rds: $(03)/bias.r $(03)/mcmc_wrong_prior.rds $(03)/mcmc_right_prior.rds 00-bias-functions.r
	$(slurmr) $(03)/bias.r
