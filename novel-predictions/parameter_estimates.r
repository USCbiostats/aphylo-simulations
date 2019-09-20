library(aphylo)
library(coda)

traceplot2 <- function(x, ylim = c(0,1), what = NULL) {
  x <- x$hist
  
  if (!is.null(what))
    x <- what(x)
  
  k <- ifelse(is.mcmc.list(x), ncol(x[[1]]), ncol(x))
  op <- par(mfrow = c(4,2), mai = c(0.3366, 0.2706, 0.2706, 0))
  traceplot(x, ylim = ylim)
  par(op)
}

# Reading the data
candidate_trees <- readRDS("data/candidate_trees_inferred.rds")

# Fitting the model
dat <- do.call(c, candidate_trees)
dat <- unlist(lapply(dat, function(d) lapply(1:Nann(d), function(i) d[i])), recursive = FALSE)
dat <- do.call(c, dat)

# Selecting using balance
b <- balance_ann(dat)
dat <- dat[(b >= quantile(b, .25))]

ans_mle <- aphylo_mle(dat ~ psi + mu_d + mu_s + Pi)

kernel_ram <- kernel_new(
  proposal = function(env) {
    
    # Making proposal
    U      <- mvtnorm::rmvnorm(1, mean = Mu, sigma = diag(7))[1,]
    theta1 <- env$theta0 + (Sigma %*% U)[,1]
    
    # Updating acceptance rate
    if (env$i > 2) {
      if (arate > env$i)
        arate <<- 0
      arate <<- arate + with(env, ans[i-1,] != ans[i-2,])
    }
    
    # Updating the variance covariance matrix
    if (!(env$i %% 50)) {

      Sigma <- Sigma %*% (
        diag(7) + min(1, (env$i)^(-2/3)*7)*(arate/env$i - .24) * U %*% t(U) /
          norm(rbind(U), "2")^2
        ) %*% t(Sigma)
      Sigma <<- chol(Sigma)
      
      # Sigma <<- cov(env$ans[1L:(env$i-1),])*2.38^2/7 + diag(7)*1e-6
      
      message("Updating the variance co-variance matrix, iteration #",  env$i)
      print(Sigma)
    }
    
    # Reflection
    reflect_on_boundaries(theta1, lb = lb, ub = ub, which = 1:7)
    
  },
  Sigma = diag(7)*1e-4,
  Mu    = rep(0, 7),
  lb    = rep(0, 7),
  ub    = rep(1, 7),
  arate = 0L
)

priors. <- function(p) 1L

mcmc. <- list(
  nchains      = 1,
  multicore    = FALSE, 
  burnin       = 0L,
  nsteps       = 10000L,
  conv_checker = NULL,
  kernel       = fmcmc::kernel_adapt(
    Sigma = chol(crossprod(vcov(ans_mle))),
    lb = 0, ub = 1, warmup = 1000, freq = 10
    ),
  thin         = 1L
)

# debug(mcmc.$kernel$proposal)
res1 <- aphylo_mcmc(
  dat ~ mu_d + mu_s + psi + Pi,
  priors = priors., control = mcmc.,
  params = coef(ans_mle)
)

res2 <- aphylo_mcmc(
  dat ~ mu_d + mu_s + psi + Pi,
  priors = priors., control = mcmc.,
  params = rep(.9, 7)
)

res3 <- aphylo_mcmc(
  dat ~ mu_d + mu_s + psi + Pi,
  priors = priors., control = mcmc.,
  params = rep(.5, 7)
)

saveRDS(res, "novel-predictions/parameter_estimates_final.rds")

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