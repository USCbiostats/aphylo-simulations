#' ---
#' title: "Convergence"
#' author: "George G Vega Yon"
#' date: "`r paste('This version:', Sys.time())`"
#' output: pdf_document
#' ---
#' 

#+ setup, echo=FALSE
knitr::opts_chunk$set(echo = FALSE)

library(sluRm)

#+ data-loading, cache=TRUE
source("simulations/00-global-parameters.r")
source("simulations/00-bias-functions.r")

dat <- readRDS("simulations/dgp.rds")

library(ggplot2)
library(magrittr)
library(dplyr)
library(aphylo)

# Computing bias ---------------------------------------------------------------
bias_MCMC_right <- bias_calc("simulations/02-missinglabel/mcmc_right_prior_estimates.rds", dat)
bias_MCMC_wrong <- bias_calc("simulations/02-missinglabel/mcmc_wrong_prior_estimates.rds", dat)

bias <- rbind(
  data.frame(Prior = "Wrong", bias_MCMC_wrong),
  data.frame(Prior = "Right", bias_MCMC_right)
)

# Categorial variables ---------------------------------------------------------

# Missings
bias$miss_tag <- interval_tags(bias$Missing, seq(0, 1, length.out = 5))

saveRDS(bias, file = "simulations/02-missinglabel/bias.rds")


