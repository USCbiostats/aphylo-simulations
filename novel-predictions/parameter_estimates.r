library(aphylo)
library(coda)

traceplot2 <- function(x) {
  x <- x$hist
  k <- ifelse(is.mcmc.list(x), ncol(x[[1]]), ncol(x))
  op <- par(mfrow = c(4,2), mai = c(0.3366, 0.2706, 0.2706, 0))
  traceplot(x, ylim = c(0,1))
  par(op)
}

# Reading the data
candidate_trees <- readRDS("data/candidate_trees.rds")

# Fitting the model
dat <- do.call(c, candidate_trees)

# dat <- dat[Ntip(dat) < 1000L]

priors. <- function(p) 1L

mcmc. <- list(
  nchains      = 1,
  multicore    = FALSE, 
  burnin       = 1000,
  nsteps       = 2000,
  conv_checker = NULL,
  kernel       = fmcmc::kernel_normal_reflective(
    scale     = .05,
    ub        = 1,
    lb        = 0,
    fixed     = FALSE,
    scheme    = "random"
  )
)
# 
# res0 <- aphylo_mcmc(
#   dat ~ mu_d + mu_s + psi + Pi,
#   priors = priors., control = mcmc.
#   )
# 
# res0b <- aphylo_mcmc(
#   dat ~ mu_d + mu_s + Pi,
#   priors = priors., control = mcmc.
# )
# 
# 
# mcmc.$burnin <- 0L
# mcmc.$nsteps <- 1e4
# 
# res1 <- aphylo_mcmc(
#   dat ~ mu_d + mu_s + psi + Pi,
#   priors = priors., control = mcmc.,
#   params = summary(res0$hist)$quantiles[,"50%"]
# )
# 
# # Adjust the scale using the acceptance rate
# arate <- 1 - rejectionRate(res1$hist)
# mcmc.$kernel$scale <- mcmc.$kernel$scale/(1 + (.44 - arate))^2
# 
# saveRDS(res1, "novel-predictions/parameter_estimates_baseline.rds")

res1 <- readRDS("novel-predictions/parameter_estimates_baseline.rds")


# Using adaptation as described by Handbook of MCMC p. 104
library(fmcmc)
kernel_mvn <- kernel_new(
  proposal = function(env) {
    
    # Updating the variance covariance matrix
    if (env$i > 1000L && !(env$i %% 500)) {
      Sigma <<- cov(env$ans[(env$i-999):(env$i-1),])*2.38^2/7
      message("Updating the variance co-variance matrix")
      print(Sigma)
    }
    
    theta1 <- env$theta0 +
      as.vector(mvtnorm::rmvnorm(1, mean = Mu, sigma = Sigma))
    reflect_on_boundaries(theta1, lb = lb, ub = ub, which = 1:7)
    
  },
  Sigma = cov(do.call(rbind, res1$hist))*2.38^2/7,
  Mu    = rep(0, 7),
  lb    = rep(0, 7),
  ub    = rep(1, 7)
)

mcmc.$kernel <- kernel_mvn
mcmc.$nsteps <- 10000
mcmc.$thin   <- 1

res2 <- aphylo_mcmc(
  dat ~ mu_d + mu_s + psi + Pi,
  priors = priors., control = mcmc.,
  params = summary(res1$hist)$quantiles[,"50%"]
)

saveRDS(res3, "novel-predictions/parameter_estimates_final.rds")

# 
# set.seed(11)
# x <- raphylo(20)
# x <- rdrop_annotations(x, .2)
# 
# psi  <- c(.1, .02)
# mu_d <- c(.7, .5)
# mu_s <- psi
# Pi   <- .5
# eta  <- c(.5, .5)
# 
# LogLike(x, psi = psi, mu_d = mu_d, mu_s = mu_s, eta = eta, Pi = Pi)$ll
# LogLike(x, psi = psi, mu_d = mu_d, mu_s = mu_s, eta = c(1,1), Pi = Pi)$ll


Updating the variance co-variance matrix
psi0          psi1         mu_d0         mu_d1
psi0   3.332543e-03 -5.778549e-07 -2.110127e-04  2.684601e-03
psi1  -5.778549e-07  3.846269e-05  2.749487e-05 -5.597484e-05
mu_d0 -2.110127e-04  2.749487e-05  3.432072e-04 -1.794077e-04
mu_d1  2.684601e-03 -5.597484e-05 -1.794077e-04  3.558301e-03
mu_s0 -1.985363e-04 -9.971792e-06  3.072210e-05  6.527238e-05
mu_s1  2.052502e-05 -7.149094e-06 -2.644419e-06  3.755578e-05
Pi    -2.330544e-03  7.977977e-05  3.759017e-04 -3.287066e-03
mu_s0         mu_s1            Pi
psi0  -1.985363e-04  2.052502e-05 -2.330544e-03
psi1  -9.971792e-06 -7.149094e-06  7.977977e-05
mu_d0  3.072210e-05 -2.644419e-06  3.759017e-04
mu_d1  6.527238e-05  3.755578e-05 -3.287066e-03
mu_s0  1.576261e-04  2.607852e-05 -2.167044e-04
mu_s1  2.607852e-05  1.723467e-05 -9.271901e-05
Pi    -2.167044e-04 -9.271901e-05  1.139612e-02