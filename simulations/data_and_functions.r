rm(list = ls())

library(aphylo)
library(coda)

# Number of simulations (samples) to draw
nsim <- 13096

# Function to draw parameters:
# This creates 5 + 1 parameters, namely
#  - psi0 Probability of mislabel a 0
#  - psi1 Probability of mislabel a 1
#  - mu0  Probability of gain a 1
#  - mu1  Probability of loss a 1
#  - Pi   Root node probability
#  - % of missings
draw_par <- function() {
  structure(
    c(rbeta(5, 1, 20), runif(1, .1, .5)),
    names = c("psi0", "psi1", "mu0", "mu1", "Pi", "missing")
    )
}

# Function to simulate annotations on a given tree
read_and_sim <- function(fn, p) {
  # Reading
  dat <- ape::read.tree(fn, nmax=1)
  
  # Simulating
  sim_annotated_tree(
    tree = dat$edge,
    psi = p[1:2],
    mu  = p[3:4],
    Pi  = p[5]
    )
}

# Function to estimate model using mle
mle_lite <- function(dat, abc, priors = NULL) {
  # Try to estimate the model
  ans <- if (abc) tryCatch(phylo_mle(dat, method="ABC", priors=priors), error = function(e) e)
    else tryCatch(phylo_mle(dat, priors = priors), error = function(e) e)
  
  
  # If it didn't worked, then return error
  if (inherits(ans, "error")) return(ans)
  
  # ans[c("par", "ll", "counts", "convergence", "message", "method", "varcovar")]
  ans
}

# Function to estimate model using MCMC
mcmc_lite <- function(dat, par, nbatch = 1e4L, nchains = 4L, burnin=1e3, thin=20, priors = NULL) {
  
  # Try to estimate the model
  ans <- tryCatch(phylo_mcmc(
    par, dat,
    control = list(nbatch = nbatch, nchains=nchains, burnin = burnin, thin=thin),
    priors = priors
    ),
    error = function(e) e
    )
  
  # If it didn't worked, then return error
  if (inherits(ans, "error")) return(ans)
  
  # ans[["dat"]] <- NULL
  
  return(ans)
  
  
}

# Function to sample a tree, and simulate annotations on it
tree_files <- list.dirs("../PANTHER11.1/books", full.names = TRUE, recursive = FALSE)
tree_files <- paste0(tree_files, "/tree.tree")

# Checking which ones exists
tree_files_exists <- sapply(tree_files, file.exists)
table(tree_files_exists) # 13096 TRUE

# Sampling files and getting the seed
set.seed(1133)
sample_files <- sample(tree_files, nsim)
parameters   <- lapply(1:nsim, function(x) draw_par())

dat          <- Map(function(fn, p) read_and_sim(fn, p), fn = sample_files, p = parameters)
dat_obs      <- Map(function(d, p) drop_data(d, p[["missing"]]), d = dat, p = parameters)

x <- lapply(dat, "[[", "annotations")
x <- lapply(x, table)
x <- lapply(x, names)
table(unlist(x))

save.image("playground/simulations/data_and_functions.rda")


