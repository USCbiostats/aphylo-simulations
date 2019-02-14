# Extract and compile Gelman tests
library(aphylo)
library(amcmc)
library(coda)

# This function extract the gelman statistics for the ith element of 
# -mcmc_right_prior_estimates_sample-
get_gelman <- function(x, obj) {
  
  # Trying to compute gelman test
  ans <- tryCatch(gelman.diag(obj[[x]][["hist"]]),
                  error = function(e) e)
  
  # Returning if error (with NA)
  if (inherits(ans, "error")) {
    print(ans)
    message("Element number -",x,"- has errors.")
    return(rep(NA, 5*2 + 1))
  }
  
  # Coercing into a vector
  dn <- rownames(ans$psrf)
  dn <- c(
    sprintf("Point est. %s", dn),
    sprintf("Upper C.I. %s", dn),
    "mpsrf"
  )
  
  # Nice named matrix
  matrix(
    c(ans$psrf[,1], ans$psrf[,2], ans$mpsrf),
    byrow = TRUE, nrow = 1,
    dimnames = list(x, dn)
  )
}


library(parallel)