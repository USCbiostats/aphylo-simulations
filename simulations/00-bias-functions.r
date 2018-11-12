# Function to calculate the coverage 
#' @param x An object of class `ahpylo_estimates`
#' @param par0 The truth vector of parameters.
#' @param pcent Type I error.
coverage <- function(x, par0, pcent) {
  quant <- do.call(rbind, x$hist)
  quant <- apply(quant, 2, quantile, c(pcent/2, 1 - pcent/2))
  
  (par0 >= quant[1, ]) & (par0 <= quant[2, ])
}

# Function to measure bias
bias_calci <- function(x, par0, tree) {
  
  if (!length(x)) return(NULL)
  
  # Names of the objects that will be stored
  # vnames <- names(coef(x))
  cnames <- c(
    "TreeSize",
    "NLeafs",
    "Missings",
    paste(vnames, "pop", sep="_"),
    paste(vnames, "estimated", sep="_"),
    paste(vnames, "bias", sep="_"),
    paste(vnames, "var", sep="_"),
    paste(vnames, "covered95", sep="_")
  )
  
  # Checking if it was able to solve it or not
  if (inherits(x, "error"))
    return(structure(rep(NA, length(cnames)), names = cnames))
  
  # Number of offspring and internal nodes
  treesize <- length(tree$offspring)
  nleafs   <- nrow(tree$tip.annotation)
  
  vrs <- diag(x$varcovar)
  if (length(vrs) == 4)
    vrs <- c(vrs, NA)
  
  structure(
    c(
      treesize,
      nleafs,
      par0["drop"],
      par0[vnames],
      coef(x),
      coef(x) - par0[vnames],
      vrs,
      coverage(x, par0[vnames], .05)
    ),
    names = cnames
  )
  
}

bias_calc <- function(fn, objname) {
  
  # Retrieving parameters and true trees
  parameters <- lapply(dat, "[[", "par")
  trees      <- lapply(dat, "[[", "atree")
  
  
  # Reading the results
  env <- new.env()
  assign(objname, readRDS(fn), envir = env)
  
  # Which ones are complete?
  ids <- which(sapply(env[[objname]], length) > 0)
  
  message("Computing bias ...", appendLF = FALSE)
  ans <- parallel::mcMap(bias_calci, env[[objname]][ids], parameters[ids], trees[ids], mc.cores=10L)
  ans <- do.call(rbind, ans)
  message("done.")
  
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