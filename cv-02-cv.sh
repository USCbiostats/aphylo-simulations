#!/bin/sh
#SBATCH --job-name="00-cv"
#SBATCH --output="hpc-logs/cv-02-slurm%A.out"
#SBATCH --time="12:00:00"
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
#SBATCH --mail-type=ALL
/usr/usc/R/3.4.4/lib64/R/bin/Rscript --vanilla \
	cv/cross-validation.r
