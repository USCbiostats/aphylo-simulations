#!/bin/sh
#SBATCH --job-name="03-misslabel-right-prior"
#SBATCH --output="03-misslabel-right-prioir-slurm%A.out"
#SBATCH --time="12:00:00"
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
#SBATCH --mail-type=ALL
/usr/usc/R/3.4.4/lib64/R/bin/Rscript --vanilla \
	simulations/03-misslabel/mcmc_right_prior_estimates.r
