#!/bin/sh
#SBATCH --job-name="01-gold-standard-wrong-prior"
#SBATCH --output="01-golds-standard-wrong-prioir-slurm%A.out"
#SBATCH --time="12:00:00"
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
/usr/usc/R/3.4.0/lib64/R/bin/Rscript --vanilla \
	simulations/01-gold-standard/mcmc_wrong_prior_estimates.r 

