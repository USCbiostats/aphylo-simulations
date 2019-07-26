library(aphylo)
library(sluRm)

source("global-paths.r")
source("simulations/00-global-parameters.r")
trees <- readRDS("data/candidate_trees.rds")

opts_sluRm$set_opts(account="lc_pdt", partition="thomas")
opts_sluRm$verbose_on()

# Unlisting functions
trees <- lapply(names(trees), function(i) {
  
  tree <- trees[[i]]
  ans  <- vector("list", Nann(tree))
  
  for (j in seq_along(ans))
    ans[[j]] <- tree[j]
  
  names(ans) <- paste0(i, "-", 1:length(ans))
  ans
  
})

trees <- unlist(trees, recursive = FALSE)

f <- function(d) {
  
  model0 <- d ~ psi + Pi + mu_d
  environment(model0) <- environment()
  
  model1 <- d ~ psi + Pi + mu_d + mu_s
  environment(model1) <- environment()
  
  ans0 <- tryCatch(aphylo_cv(
    model0,
    control = list(
      nsteps       = 1e5L,
      nchains      = 4L,
      burnin       = 2e4L,
      thin         = 50L,
      multicore    = FALSE,
      conv_checker = amcmc::gelman_convergence(5e3, 1.05)
    ),
    priors = bprior(2,20)
    ), error = function(e) e)
  
  ans1 <- tryCatch(aphylo_cv(
    model1,
    control = list(
      nsteps       = 1e5L,
      nchains      = 4L,
      burnin       = 2e4L,
      thin         = 50L,
      multicore    = FALSE,
      conv_checker = amcmc::gelman_convergence(5e3L, 1.05)
    ),
    priors = function(p) {
      dbeta(
        x      = p[c("mu_d0", "mu_d1", "mu_s0", "mu_s1", "psi0", "psi1", "Pi")],
        shape1 = c(20, 20, 2, 2, 2, 2, 2),
        shape2 = c(2, 2, 20, 20, 20, 20, 20)
        )
      }
  ), error = function(e) e)
  
  list(ans0 = ans0, ans1 = ans1)
  
}


ans_slurm <- Slurm_lapply(
  trees, f, njobs=50, mc.cores=4L,
  job_name="aphylo-cross-validation",
  plan = "wait"
  )
saveRDS(ans_slurm, "cv/cross_validation-job.rds")
saveRDS(Slurm_collect(ans_slurm), "cv/cross_validation.rds")


