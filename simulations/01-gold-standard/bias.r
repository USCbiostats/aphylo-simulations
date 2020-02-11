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
source("00-global-parameters.r")
source("00-bias-functions.r")

dat <- readRDS("dgp.rds")

library(ggplot2)
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


