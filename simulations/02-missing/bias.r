#' ---
#' title: "Convergence"
#' author: "George G Vega Yon"
#' date: "`r paste('This version:', Sys.time())`"
#' output: pdf_document
#' ---
#' 

#+ setup, echo=FALSE
knitr::opts_chunk$set(echo = FALSE)

#+ data-loading, cache=TRUE
dat <- readRDS("simulations/dgp.rds")


library(ggplot2)
library(magrittr)

coverage <- function(x, par0, pcent) {
  quant <- do.call(rbind, x$hist)
  quant <- apply(quant, 2, quantile, c(pcent/2, 1 - pcent/2))
  
  (par0 >= quant[1, ]) & (par0 <= quant[2, ])
}

# Function to measure bias
bias_calci <- function(x, par0, tree) {

  if (!length(x)) return(NULL)

  # Names of the objects that will be stored
  parnames <- c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1", "Pi") %>%
    intersect(names(coef(x)))
  cnames <- c(
    "TreeSize",
    "NLeafs",
    "Missings",
    parnames,
    paste(parnames, "bias", sep="_"),
    paste(parnames, "var", sep="_"),
    paste(parnames, "covered95", sep="_")
  )
  
  # Checking if it was able to solve it or not
  if (inherits(x, "error"))
    return(structure(rep(NA, length(cnames)), names = cnames))

  # Number of offspring and internal nodes
  treesize <- length(tree$offspring)
  nleafs   <- nrow(tree$tip.annotation)

  vrs <- diag(vcov(x))
  if (length(vrs) == 4)
    vrs <- c(vrs, NA)
  
  structure(
    c(
      treesize,
      nleafs,
      par0["drop"],
      x$par,
      x$par - par0[parnames],
      vrs,
      coverage(x, par0[parnames], .05)
      ),
    names = cnames
    )

}

bias_calc <- function(x, parameters, trees) {
  
  ids <- which(sapply(x, length) > 0)

  ans <- parallel::mcmapply(
    bias_calci, x[ids], parameters[ids], trees, mc.cores = 9
    )
  
  ans <- do.call(rbind, ans)

  ans <- cbind(index = ids, ans)

  ans[complete.cases(ans[,-ncol(ans)]),,drop=FALSE]
}

# Creates nice interval tags in the form of [a,b)...[a,z] (last one closed).
# All numbers must be within the interval
interval_tags <- function(x, marks) {

  # Generating labels
  n <- length(marks)
  l <- c(sprintf("[%.1f, %.1f)", marks[-n][-(n-1)], marks[-n][-1]),
    sprintf("[%.1f, %.1f]", marks[n-1], marks[n])
  )

  # Finding intervals
  x <- findInterval(x, marks, rightmost.closed = TRUE)
  factor(x, levels = 1:length(l), labels = l)

}


# Computing bias ---------------------------------------------------------------
# bias_MLE <- bias_calc("simulations/mle_estimates.rda", "ans_MLE")
# bias_MLE2 <- bias_calc("simulations/mle_estimates2.rda", "ans_MLE")
# bias_MAP <- bias_calc("simulations/map_estimates.rda", "ans_MAP")
# bias_MAP_wrong <- bias_calc("simulations/map_wrong_prior_estimates.rda", "ans_MAP_wrong_prior")
mcmc_right_prior_estimates <- readRDS("simulations/02-missing/mcmc_right_prior_estimates.rds")
bias_MCMC_right <- #bias_calc("simulations/02-missing/mcmc_right_prior_estimates.rda", "ans_MCMC_right_prior")
  bias_calc(mcmc_right_prior_estimates, lapply(dat, "[[", "par"), lapply(dat, "[[", "atree"))
# bias_MCMC_wrong <- bias_calc("simulations/02-missing/mcmc_wrong_prior_estimates.rda", "ans_MCMC_wrong_prior")

# Checking solved solutions
# common_solutions <- bias_MLE[,"index"]
# common_solutions <- intersect(common_solutions, bias_MAP_wrong[,"index"])
# common_solutions <- intersect(common_solutions, bias_MCMC_right[,"index"])
# common_solutions <- intersect(common_solutions, bias_MCMC_wrong[,"index"])

bias <- rbind(
  # data.frame(Method = "MLE", bias_MLE),
  # data.frame(Method = "MLE2", bias_MLE2)
  # data.frame(Method = "MAP", bias_MAP),
  # data.frame(Method = "MAP wrong", bias_MAP_wrong),
  # data.frame(Prior = "Wrong", bias_MCMC_wrong),
  data.frame(Prior = "Right", bias_MCMC_right)
)

# Categorial variables ---------------------------------------------------------

# Missings
bias$miss_tag <- interval_tags(bias$Missings, 1:5/ 10)

# Tree size
bias$size_tag <- interval_tags(bias$TreeSize, quantile(bias$TreeSize, na.rm = TRUE))

# NLeafs/TreeSize
bias$PropLeafs <- with(bias, NLeafs/TreeSize)
bias$PropLeafs_tag <- interval_tags(bias$PropLeafs, quantile(bias$PropLeafs, na.rm=TRUE))

save(bias, file = "simulations/02-missing/bias.rda")


