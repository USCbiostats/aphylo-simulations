#!/bin/sh
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=03-bias
#SBATCH --mail-type=ALL
#SBATCH --mail-user=g.vegayon@gmail.com

#+ setup, echo=FALSE
knitr::opts_chunk$set(echo = FALSE)

library(slurmR)

#+ data-loading, cache=TRUE
source("00-global-parameters.r")
source("00-bias-functions.r")

dat <- readRDS("dgp.rds")

library(magrittr)
library(dplyr)
library(aphylo)

# Computing bias ---------------------------------------------------------------
parnames <- c("psi0", "psi1", "mu_d0", "mu_d1", "mu_s0", "mu_s1", "Pi")
bias_MCMC_right <- bias_calc("03-misslabel/mcmc_right_prior.rds", dat)
bias_MCMC_wrong <- bias_calc("03-misslabel/mcmc_wrong_prior.rds", dat)

bias <- rbind(
  data.frame(Prior = "Wrong", bias_MCMC_wrong),
  data.frame(Prior = "Right", bias_MCMC_right)
)

# Categorial variables ---------------------------------------------------------

# Missings
bias$miss_tag <- interval_tags(bias$Missing, seq(0, 1, length.out = 5))

saveRDS(bias, file = "03-misslabel/bias.rds")


