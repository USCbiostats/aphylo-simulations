slurmr := ~/./slurmr
.PHONY: all
01:=01-gold-standard
02:=02-missinglabel
03:=03-misslabel

all: */mcmc_right_prior.rds */mcmc_wrong_prior.rds 

# Data generating process
dgp.rds: dgp.r
	$(slurmr) dgp.r

# Running the gold standard ----------------------------------------------------
$(01)/mcmc_right_prior.rds: $(01)/mcmc_right_prior.r dgp.rds
	$(slurmr) $(01)/mcmc_right_prior.r
$(01)/mcmc_wrong_prior.rds: $(01)/mcmc_wrong_prior.r dpg.rds
	$(slurmr) $(01)/mcmc_right_prior.r

# Running simulations with missing data ----------------------------------------
$(02)/mcmc_right_prior.rds: $(02)/mcmc_right_prior.r dgp.rds
	$(slurmr) $(02)/mcmc_right_prior.r
$(02)/mcmc_wrong_prior.rds: $(02)/mcmc_wrong_prior.r dpg.rds
	$(slurmr) $(02)/mcmc_right_prior.r

# Running simulations with missing data and mis labels -------------------------
$(03)/mcmc_right_prior.rds: $(03)/mcmc_right_prior.r dgp.rds
	$(slurmr) $(03)/mcmc_right_prior.r
$(03)/mcmc_wrong_prior.rds: $(03)/mcmc_wrong_prior.r dpg.rds
	$(slurmr) $(03)/mcmc_right_prior.r

