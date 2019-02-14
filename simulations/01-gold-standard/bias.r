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
vnames <- c("mu0", "mu1", "Pi")
bias_MCMC_right <- bias_calc("simulations/01-gold-standard/mcmc_right_prior_estimates.rds", "ans_MCMC_right_prior")
bias_MCMC_wrong <- bias_calc("simulations/01-gold-standard/mcmc_wrong_prior_estimates.rds", "ans_MCMC_wrong_prior")

bias <- rbind(
  data.frame(Prior = "Wrong", bias_MCMC_wrong),
  data.frame(Prior = "Right", bias_MCMC_right)
)

# Categorial variables ---------------------------------------------------------

# Missings
bias$Missing  <- 0L # This case has no missigness
bias$miss_tag <- "No missings"

# Tree size
bias$size_tag <- interval_tags(bias$TreeSize, quantile(bias$TreeSize, na.rm = TRUE))

# NLeafs/TreeSize
bias$PropLeafs <- with(bias, NLeafs/TreeSize)
bias$PropLeafs_tag <- interval_tags(bias$PropLeafs, quantile(bias$PropLeafs, na.rm=TRUE))

saveRDS(bias, file = "simulations/01-gold-standard/bias.rds")


