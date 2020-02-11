#!/bin/sh
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=01-bias
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
bias_MCMC_right <- bias_calc("01-gold-standard/mcmc_right_prior.rds", dat)
bias_MCMC_wrong <- bias_calc("01-gold-standard/mcmc_wrong_prior.rds", dat)

bias <- rbind(
  data.frame(Prior = "Wrong", bias_MCMC_wrong),
  data.frame(Prior = "Right", bias_MCMC_right)
)

# Categorial variables ---------------------------------------------------------

# Missings
bias$Missing  <- 0L # This case has no missigness
bias$miss_tag <- "No missings"

saveRDS(bias, file = "01-gold-standard/bias.rds")


