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

res0 <- aphylo_mcmc(
  dat ~ mu_d + mu_s + psi + Pi,
  priors = priors., control = mcmc.
  )


mcmc.$burnin <- 0L
mcmc.$nsteps <- 1e4

res1 <- aphylo_mcmc(
  dat ~ mu_d + mu_s + psi + Pi,
  priors = priors., control = mcmc.,
  params = summary(res0$hist)$quantiles[,"50%"]
)

# Adjust the scale using the acceptance rate
arate <- 1 - rejectionRate(res1$hist)
mcmc.$kernel$scale <- mcmc.$kernel$scale/(1 + (.44 - arate))^2

saveRDS(res1, "novel-predictions/parameter_estimates_baseline.rds")

# res2 <- readRDS("novel-predictions/parameter_estimates.rds")


# Using adaptation as described by Handbook of MCMC p. 104
library(fmcmc)
kernel_mvn <- kernel_new(
  proposal = function(env) {
    theta1 <- env$theta0 +
      as.vector(mvtnorm::rmvnorm(1, mean = Mu, sigma = Sigma))
    reflect_on_boundaries(theta1, lb = lb, ub = ub, which = 1:7)
    
  },
  Sigma = cov(do.call(rbind, res2$hist))*2.38^2/7,
  Mu    = rep(0, 7),
  lb    = rep(0, 7),
  ub    = rep(1, 7)
)

mcmc.$kernel <- kernel_mvn
mcmc.$nsteps <- 10000
mcmc.$thin   <- 1

res3 <- aphylo_mcmc(
  dat ~ mu_d + mu_s + psi + Pi,
  priors = priors., control = mcmc.,
  params = summary(res2$hist)$quantiles[,"50%"]
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
