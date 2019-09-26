library(aphylo)
library(coda)

shrink_towards_half <- function(x, margin=.01) {
  
  x[x < (.5 - margin)] <- x[x < (.5 - margin)] + margin
  x[x > (.5 + margin)] <- x[x > (.5 + margin)] - margin
  
  x
}

# Common parameters
prior.  <- bprior(c(2,2,9,9,2,2,2,2,2), c(9,9,2,2,9,9,9,9,9))
warmup. <- 1000
freq.   <- 10
lb.     <- 1e-5
ub.     <- 1 - 1e-5

mcmc. <- list(
  nchains      = 4L,
  multicore    = FALSE, 
  burnin       = 500L,
  nsteps       = 10000L,
  conv_checker = NULL,
  kernel       = fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup., freq = freq.),
  thin         = 10L
)

set.seed(1362)

# Case 1: (Almost) Fully annotated trees ---------------------------------------

# Data preprocessing
fully_annotated <- readRDS("data/candidate_trees_inferred.rds")
fully_annotated <- do.call(c, fully_annotated)
fully_annotated <- unlist(lapply(fully_annotated, function(d) {
  lapply(1:Nann(d), function(i) d[i])
  }), recursive = FALSE)
fully_annotated <- do.call(c, fully_annotated)

# Selecting using balance
b <- balance_ann(fully_annotated)
fully_annotated <- fully_annotated[(b >= quantile(b, .2))]

# 1.A: No prior
ans_mle_fully_annotated_no_prior <- aphylo_mle(
  fully_annotated ~ psi + mu_d + mu_s + Pi
)

message("Fully annotated: MLE No prior done.")

saveRDS(
  ans_mle_fully_annotated_no_prior,
  "novel-predictions/parameter_estimates_mle_fully_annotated_no_prior.rds"
  )

ans_mcmc_fully_annotated_no_prior <- aphylo_mcmc(
  fully_annotated ~ psi + mu_d + mu_s + Pi,
  params  = shrink_towards_half(coef(ans_mle_fully_annotated_no_prior)),
  control = mcmc.
)

message("Fully annotated: MCMC No prior done.")

saveRDS(
  ans_mcmc_fully_annotated_no_prior,
  "novel-predictions/parameter_estimates_mcmc_fully_annotated_no_prior.rds"
  )
# 1.B: Prior
ans_mle_fully_annotated_prior <- aphylo_mle(
  fully_annotated ~ psi + mu_d + mu_s + Pi, priors = prior.
)

message("Fully annotated: MLE prior done.")

saveRDS(
  ans_mle_fully_annotated_prior,
  "novel-predictions/parameter_estimates_mle_fully_annotated_prior.rds"
  )
mcmc.$kernel <- fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup., freq = freq.)
ans_mcmc_fully_annotated_prior <- aphylo_mcmc(
  fully_annotated ~ psi + mu_d + mu_s + Pi,
  priors = prior.,
  params  = shrink_towards_half(coef(ans_mle_fully_annotated_prior)),
  control = mcmc.
)

message("Fully annotated: MCMC prior done.")

saveRDS(
  ans_mcmc_fully_annotated_prior,
  "novel-predictions/parameter_estimates_mcmc_fully_annotated_prior.rds"
  )
# Case 2: (Almost) Fully annotated trees ---------------------------------------

partially_annotated <- readRDS("data/candidate_trees.rds")
partially_annotated <- do.call(c, partially_annotated)
partially_annotated <- unlist(lapply(partially_annotated, function(d) {
  lapply(1:Nann(d), function(i) d[i])
}), recursive = FALSE)
partially_annotated <- do.call(c, partially_annotated)

# 1.A: No prior
ans_mle_partially_annotated_no_prior <- aphylo_mle(
  partially_annotated ~ psi + mu_d + mu_s + Pi
)

message("Partially annotated: MLE No prior done.")

saveRDS(
  ans_mle_partially_annotated_no_prior,
  "novel-predictions/parameter_estimates_mle_partially_annotated_no_prior.rds"
  )
mcmc.$kernel <- fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup., freq = freq.)
ans_mcmc_partially_annotated_no_prior <- aphylo_mcmc(
  partially_annotated ~ psi + mu_d + mu_s + Pi,
  params  = shrink_towards_half(coef(ans_mle_partially_annotated_no_prior)),
  control = mcmc.
)

message("Partially annotated: MCMC No prior done.")

saveRDS(
  ans_mcmc_partially_annotated_no_prior,
  "novel-predictions/parameter_estimates_mcmc_partially_annotated_no_prior.rds"
  )
# 1.B: Prior
ans_mle_partially_annotated_prior <- aphylo_mle(
  partially_annotated ~ psi + mu_d + mu_s + Pi, priors = prior.
)

message("Partially annotated: MLE prior done.")

saveRDS(
  ans_mle_partially_annotated_prior,
  "novel-predictions/parameter_estimates_mle_partially_annotated_prior.rds"
  )
mcmc.$kernel <- fmcmc::kernel_adapt(lb = lb., ub = ub., warmup = warmup., freq = freq.)
ans_mcmc_partially_annotated_prior <- aphylo_mcmc(
  partially_annotated ~ psi + mu_d + mu_s + Pi,
  priors = prior.,
  params  = shrink_towards_half(coef(ans_mle_partially_annotated_prior)),
  control = mcmc.
)

message("Partially annotated: MCMC prior done.")

saveRDS(
  ans_mcmc_partially_annotated_prior,
  "novel-predictions/parameter_estimates_mcmc_partially_annotated_prior.rds"
  )
