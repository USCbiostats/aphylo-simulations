#!/bin/sh
#SBATCH --job-name="03-pub-bias-right-prior"
#SBATCH --output="03-pub-bias-right-prior-slurm%A.out"
#SBATCH --time="5:00:00"
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
/usr/usc/R/3.4.4/lib64/R/bin/Rscript --vanilla \
  simulations/03-pub-bias/mcmc_right_prior_estimates.r 
