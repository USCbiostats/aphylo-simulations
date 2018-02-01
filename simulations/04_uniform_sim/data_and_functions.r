rm(list = ls())

library(aphylo)
library(coda)

# Number of simulations (samples) to draw
nsim   <- 13096

mcmc.nbatch  <- 2e4
mcmc.burnin  <- 5e3
mcmc.thin    <- 20
mcmc.nchains <- 4


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
    c(rbeta(2, 1, 39), rbeta(3, 1, 19), runif(1, .1, .5)),
    names = c("psi0", "psi1", "mu0", "mu1", "Pi", "missing")
    )
}

# Function to simulate annotations on a given tree
read_and_sim <- function(fn, p) {
  # Reading
  dat <- ape::read.tree(fn, nmax=1)
  
  # Simulating
  ans <- try(sim_annotated_tree(
    tree = dat$edge,
    psi = p[1:2],
    mu  = p[3:4],
    Pi  = p[5],
    informative = FALSE
    ))
  
  
}

# Function to estimate model using mle
mle_lite <- function(par, dat, abc, priors = NULL) {
  # Try to estimate the model
  ans <- if (abc) tryCatch(aphylo_mle(dat, method="ABC", priors=priors, par = par), error = function(e) e)
    else tryCatch(aphylo_mle(dat, priors = priors, par = par), error = function(e) e)
  
  
  # If it didn't worked, then return error
  if (inherits(ans, "error")) return(ans)
  
  # ans[c("par", "ll", "counts", "convergence", "message", "method", "varcovar")]
  ans
}

# Function to estimate model using MCMC
mcmc_lite <- function(
  dat,
  par,
  priors  = NULL,
  nbatch  = mcmc.nbatch,
  nchains = mcmc.nchains,
  burnin  = mcmc.burnin,
  thin    = mcmc.thin
  ) {
  
  # Try to estimate the model
  ans <- tryCatch(aphylo_mcmc(
    par, dat,
    control = list(nbatch = nbatch, nchains=nchains, burnin = burnin, thin=thin),
    priors = priors,
    check.informative = FALSE
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

dat <- lapply(1:nsim, function(x) {
  
  while (TRUE) {
    
    # Drawing parameters
    par <- draw_par()
    
    # Sampling tree
    tree <- sample(tree_files, 1)
    
    # Simulating tree and annotations
    ann_tree <- read_and_sim(tree, par)
  
    # Dropping data
    ann_tree_obs <- rdrop_annotations(ann_tree, par[6], informative = FALSE)
    
    # Is it informative?
    info_obs <- aphylo:::fast_table_using_labels(ann_tree_obs$tip.annotation, c(0, 1, 9))
    if (!prod(info_obs))
      next
    else 
      break
  }
  
  # Dropping internal node annotations
  ann_tree_obs$node.annotation[] <- 9L
  
  list(
    dat      = ann_tree,
    dat_obs = ann_tree_obs,
    par     = par,
    counts    = data.frame(
      dat     = aphylo:::fast_table_using_labels(ann_tree$tip.annotation, c(0, 1, 9)),
      dat_obs = info_obs
    )
  )
})


dat_obs    <- lapply(dat, "[[", "dat_obs")
parameters <- lapply(dat, "[[", "par")
dat        <- lapply(dat, "[[", "dat")


# How many are informative?
info <- lapply(
  lapply(dat_obs, "[[", "tip.annotation"),
  aphylo:::fast_table_using_labels, ids = c(0L, 1L, 9L)
  )

info <- apply(do.call(cbind, info), 2, prod)
stopifnot(all(info > 0))
rm(info)


x <- lapply(dat_obs, "[[", "tip.annotation")
x <- lapply(x, aphylo:::fast_table_using_labels, ids = c(0L, 1L, 9L))
x <- t(do.call(cbind, x))
x <- x/rowSums(x)
boxplot(x)

rm(x)
save.image("simulations/04_uniform_sim/data_and_functions.rda", compress = FALSE)

