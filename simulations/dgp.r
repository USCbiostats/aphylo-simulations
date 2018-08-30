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

# Path to the current panther dataset
PANTHER_PATH <- "/auto/pmd-02/pdt/pdthomas/panther/famlib/rel/PANTHER13.1_altVersion/hmmscoring/PANTHER13.1/books"
PROJECT_PATH <- "/home/rcf-proj2/pdt/vegayon/aphylo-simulations"

NSAMPLES     <- 1e4

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
sim_annotations_panther <- function(tree, par) {
  
  # Step 0: Model parameters
  par <- rbeta(
    7,
    shape1 = c( 2, 2, 2, 2, 7,18, 2),
    shape2 = c(18,18,18,18, 3, 2,18)
    )
  names(par)  <- c("psi0", "psi1", "mu0", "mu1", "eta0", "eta1", "Pi")
  par["drop"] <- runif(1, .1, .9)
  
  # 1. Simulate annotations
  atree <- sim_fun_on_tree(
    tree,
    psi = c(0, 0),
    mu  = par[c("mu0", "mu1")],
    eta = c(1, 1),
    Pi  = par["Pi"]
    )
  
  nleaf <- ape::Ntip(tree)
  atree <- new_aphylo(
    tip.annotation  = atree[1:nleaf,,drop=FALSE],
    node.annotation = atree[(nleaf+1):nrow(atree),,drop=FALSE],
    tree            = tree
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

# Generating random data w/ trees of size 50
set.seed(1)
SAMPLED_TREE_FILES <- sample(PANTHER_FILES, NSAMPLES, TRUE)
PANTHER_TREES      <- lapply(SAMPLED_TREE_FILES, read_panther)
PANTHER_TREES      <- lapply(PANTHER_TREES, "[[", "tree")

dat <- vector("list", NSAMPLES)
for (i in seq_along(PANTHER_TREES))
  dat[[i]] <- sim_annotations_panther(PANTHER_TREES[[i]])
  
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
