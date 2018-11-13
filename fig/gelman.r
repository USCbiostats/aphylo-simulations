rm(list = ls())

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

# Checkingout Right prior estimates --------------------------------------------
ans_MCMC_right_prior <- readRDS("simulations/01-gold-standard/mcmc_right_prior_estimates.rds")


# Extracting the data
N <- length(ans_MCMC_right_prior)
gelmans <- do.call(rbind, parallel::mclapply(1:N, get_gelman, obj = ans_MCMC_right_prior, mc.cores=4))



# Saving the data
save(gelmans, file = "fig/gelmans_right_prior.rda")

graphics.off()
pdf("fig/gelmans_right_prior.pdf")
oldpar <- par(no.readonly = TRUE)
par(mar = par("mar")*c(2, 1, 1, 1))
boxplot(gelmans, las = 2)#, main = "Distribution of\nGelman Statistics (right prior)")
par(oldpar)
dev.off()


# Checkingout Right prior estimates --------------------------------------------
load("simulations/mcmc_wrong_prior_estimates.rda")
cl <- makeForkCluster(8)

# Extracting the data
N <- length(ans_MCMC_right_prior)
gelmans <- do.call(rbind, parLapply(1:N, get_gelman, obj = ans_MCMC_wrong_prior))
stopCluster(cl)

# Saving the data
save(gelmans, file = "fig/gelmans_wrong_prior.rda")

graphics.off()
pdf("fig/gelmans_wrong_prior.pdf")
oldpar <- par(no.readonly = TRUE)
par(mar = par("mar")*c(2, 1, 1, 1))
boxplot(gelmans, las = 2)#, main = "Distribution of\nGelman Statistics (wrong prior)")
par(oldpar)
dev.off()

