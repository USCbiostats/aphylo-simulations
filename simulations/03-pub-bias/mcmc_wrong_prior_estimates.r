rm(list = ls())

library(aphylo)
library(sluRm)

source("simulations/00-global-parameters.r")
dat0 <- readRDS("simulations/dgp.rds") #[1:NSAMPLES]
dat0 <- lapply(dat0, "[[", "atree_o")

# Setting the seed
set.seed(111222)

mcmc.par   <- matrix(runif(5*mcmc.nchains), ncol=5)
mcmc.prior <- function(p) {
  dbeta(p, c(2, 2, 2, 2, 2), c(18, 18, 18, 18, 18))
}

job <- Slurm_lapply(
  dat0,
  mcmc_lite,
  model      = dat ~ mu + psi + Pi,
  params     = mcmc.par,
  priors     = mcmc.prior,
  nsteps     = mcmc.nsteps,
  nchains    = mcmc.nchains,
  burnin     = mcmc.burnin,
  thin       = mcmc.thin,
  njobs      = 55L,
  mc.cores   = 4L,
  multicore  = mcmc.multicore, # TRUE,
  job_name   = "03-pub-bias-wrong-prior",
  job_path   = STAGING_PATH,
  submit     = TRUE
 )

saveRDS(job, paste0(PROJECT_PATH, "/simulations/03-pub-bias/wrong-prior-estimates-job.rds"))

saveRDS(
  res <- Slurm_collect(job),
  paste0(
    PROJECT_PATH,
    "/simulations/03-pub-bias/mcmc_wrong_prior_estimates.rds"
    )
  )

# res
# job

source("simulations/summarize_predictions.r")

N <- length(ans_MCMC_wrong_prior)

# Prediction scores
pred_scores <- summarize_predictions(
  x    = ans_MCMC_wrong_prior,
  dat. = lapply(dat[1:N], "[[", "atree")
)

saveRDS(pred_scores, file = "simulations/03-pub-bias/mcmc_wrong_prior_prediction.rds", compress = FALSE)


ans <- lapply(pred_scores, function(x) {
  tryCatch(cbind(obs = x[["obs"]], rand = x[["random"]]) / x[["worse"]], error = function(e) NULL)
})
ans <- do.call(rbind, ans)

dimnames(ans) <- list(NULL, c("Model Predictions", "Random Predictions"))


graphics.off()
pdf("simulations/03-pub-bias/mcmc_wrong_prior_prediction.pdf")
# png("simulations/mcmc_wrong_prior_prediction.png")
boxplot(ans, #main = "Distribution Relative\nPrediction Scores",
        ylab = "Relative Prediction Score (0 is perfect prediction)",
        sub  = sprintf(
          "Random annotations based on a Bernoulli(%.2f)",
          attr(pred_scores, "prop_of_1s")),
        ylim = c(0,1),outline = FALSE, notch=FALSE
)
dev.off()

