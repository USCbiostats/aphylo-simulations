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

# Computing bias ---------------------------------------------------------------
vnames <- c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1","Pi")
bias_MCMC_right <- bias_calc("simulations/04-full-model/mcmc_right_prior_estimates.rds", "ans_MCMC_right_prior")
bias_MCMC_wrong <- bias_calc("simulations/04-full-model/mcmc_wrong_prior_estimates.rds", "ans_MCMC_wrong_prior")

bias <- rbind(
  data.frame(Prior = "Wrong", bias_MCMC_wrong),
  data.frame(Prior = "Right", bias_MCMC_right)
)

# Categorial variables ---------------------------------------------------------

# Missings
bias$miss_tag <- interval_tags(bias$Missing, seq(0.1, 0.9, length.out = 5))

# Tree size
bias$size_tag <- interval_tags(bias$TreeSize, quantile(bias$TreeSize, na.rm = TRUE))

# NLeafs/TreeSize
bias$PropLeafs <- with(bias, NLeafs/TreeSize)
bias$PropLeafs_tag <- interval_tags(bias$PropLeafs, quantile(bias$PropLeafs, na.rm=TRUE))

saveRDS(bias, file = "simulations/04-full-model/bias.rds")


