#!/bin/sh
#SBATCH --partition=thomas
#SBATCH --account=lc_pdt
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=01-dgp-aphylo-simulations

#' This R script contains the set of functions used for the Data Generating
#' Process for the simulation of annotated trees.
#' The following functions are included:
#' 
#' - `apply_mislabeling`: Takes a matrix with annotations and swaps 0/1 according
#'   to mislabeling probabilities.
#'   
#' - `apply_pub_bias`: Also takes a matrix with annotations, but sets some equal
#'   to 9 (`NA`) depending on the parameter `eta`.
#'   
#' - `sim_annotations_panther`: A function that takes the two previous ones,
#'   and simulates annotations (see paper).
#'   
#'   
#' The last section of the script does some testing to make sure that the functions
#' are doing what are supposed to do.

library(aphylo)
library(slurmR)

# Path to the current panther dataset
source("00-global-parameters.r")

# Data Generating Process Functions --------------------------------------------

# Mislabeling
apply_mislabeling <- function(x, psi) {
  
  n <- nrow(x)
  x[] <- ifelse(x == 9L, 9L, ifelse(psi[x + 1] < runif(n), x, 1-x))
  
  x
}

# Publication bias
apply_pub_bias <- function(x, eta) {
  
  n <- nrow(x)
  x[] <- ifelse(x == 9L, 9L, ifelse(eta[x + 1] > runif(n), x, 9L))
  
  x
}


# Core simulation
sim_annotations_panther <- function(dat, par) {
  
  tree <- dat$tree
  
  # Step 0: Model parameters
  par <- rbeta(
    9,
    shape1 = ALPHA_PAR,
    shape2 = BETA_PAR 
    )
  names(par)  <- c("psi0", "psi1", "mu_d0", "mu_d1", "mu_s0", "mu_s1", "eta0", "eta1", "Pi")
  par["drop"] <- runif(1, .1, .9)
  
  # 1. Simulate annotations
  atree <- sim_fun_on_tree(
    tree,
    # 0: duplication event, 1: other
    node.type = 1 - dat$internal_nodes_annotations$duplication,
    psi       = c(0, 0),
    mu_d      = par[c("mu_d0", "mu_d1")],
    mu_s      = par[c("mu_s0", "mu_s1")],
    eta       = c(1, 1),
    Pi        = par["Pi"]
    )
  
  nleaf <- ape::Ntip(tree)
  atree <- new_aphylo(
    tip.annotation  = atree[1:nleaf,,drop=FALSE],
    node.annotation = atree[(nleaf+1):nrow(atree),,drop=FALSE],
    tree            = tree,
    node.type       = 1 - dat$internal_nodes_annotations$duplication
    )
  
  # 2. Missing annotations
  atree_i <- atree
  atree_i$node.annotation[] <- 9L
  atree_i <- rdrop_annotations(atree_i, par["drop"], informative = FALSE)
  
  # Swapping annotaitons
  atree_m <- atree_i
  atree_m$tip.annotation <- apply_mislabeling(atree_m$tip.annotation, par[c("psi0", "psi1")])
  
  atree_o <- atree_m
  atree_o$tip.annotation <- apply_pub_bias(atree_o$tip.annotation, par[c("eta0", "eta1")])

  # Returning
  list(
    atree   = atree,
    atree_i = atree_i,
    atree_m = atree_m,
    atree_o = atree_o,
    par     = par
  )
  
}

# Evaluation of the DGP --------------------------------------------------------

# Listing panther trees
PANTHER_FILES <- list.files(PANTHER_PATH, recursive=FALSE, full.names=TRUE)
PANTHER_FILES <- paste0(PANTHER_FILES, "/tree.tree")

# REading the trees ------------------------------------------------------------
set.seed(1)
PANTHER_TREES <- Slurm_lapply(
  PANTHER_FILES,
  read_panther,
  job_name = "reading_trees",
  njobs    = 20,
  mc.cores = 4L,
  plan     = "collect"
  )

# PANTHER_TREES <- lapply(PANTHER_TREES, "[[", "tree")

dat <- vector("list", length(PANTHER_TREES)) # NSAMPLES)
for (i in seq_along(PANTHER_TREES))
  dat[[i]] <- tryCatch(
    sim_annotations_panther(PANTHER_TREES[[i]]),
    error = function(e) e)
  
# Saving
saveRDS(
  dat,
  sprintf("%s/simulations/dgp.rds", PROJECT_PATH),
  compress=FALSE
)

# Retrieving parameters
parameters <- do.call(rbind, lapply(dat, "[[", "par"))
op <- par(mfrow = c(4, 2), mar = c(3, 3, 1, 1))
for (p in colnames(parameters)) {
  
  # Histogram of the parameters
  hist(parameters[,p], main = p, breaks=50, border="transparent",
       col = "steelblue")
  
  # Drawing a line on the mean and annotating
  m <- mean(parameters[,p])
  abline(v=m, lwd=2, col="red", lty=2)
  text(x = m, y = 0, labels = sprintf("%.2f", m), pos = 4)
}
par(op)


# Comparing probabilities of publication bias

aa <- lapply(dat, function(d) {
  
  n1  <- sum(d$atree_m$tip.annotation == 1)
  n0  <- sum(d$atree_m$tip.annotation == 0)
  n19 <- sum(d$atree_o$tip.annotation == 9 & d$atree_m$tip.annotation == 1)
  n09 <- sum(d$atree_o$tip.annotation == 9 & d$atree_m$tip.annotation == 0)
  
  c("0" = n09/(n0 + 1e-10), "1" = n19/(n1 + 1e-10))
  
})

aa <- 1 - do.call(rbind, aa)
boxplot(aa, at = c(1, 2), main = "Publication bias (eta)")
text(
  x=c(1, 2), y = colMeans(aa), labels = sprintf("%.2f", colMeans(aa)))


# Comparing probabilities of mislabeling

aa <- lapply(dat, function(d) {
  
  n1  <- sum(d$atree$tip.annotation == 1)
  n0  <- sum(d$atree$tip.annotation == 0)
  n19 <- sum(d$atree_m$tip.annotation == 0 & d$atree$tip.annotation == 1)
  n09 <- sum(d$atree_m$tip.annotation == 1 & d$atree$tip.annotation == 0)
  
  c("0" = n09/(n0 + 1e-10), "1" = n19/(n1 + 1e-10))
  
})

aa <- do.call(rbind, aa)
boxplot(aa, at = c(1, 2), main = "Mislabeling")
text(
  x=c(1, 2), y = colMeans(aa), labels = sprintf("%.2f", colMeans(aa)))


