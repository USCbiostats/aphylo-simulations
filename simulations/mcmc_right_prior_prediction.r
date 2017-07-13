# rm(list = ls())

library(aphylo)

load("playground/simulations/data_and_functions.rda")
load("playground/simulations/mcmc_right_prior_estimates_sample.rda")

N <- 1000

library(parallel)
cl <- makeForkCluster(10)

# Prediction scores
pred_scores <- parLapply(cl, 1:N , function(i) { 
  try(prediction_score(mcmc_right_prior_estimates_sample[[i]], expected = dat[[i]]$annotations))
})

stopCluster(cl)
save(pred_scores, file = "playground/simulations/mcmc_right_prior_prediction.rda")

ans <- lapply(pred_scores, function(x) {
  tryCatch(cbind(obs = x[["obs"]], rand = x[["random"]]["mean"]) / x[["worse"]], error = function(e) NULL)
  })
ans <- do.call(rbind, ans)

dimnames(ans) <- list(NULL, c("Observed", "Random"))

graphics.off()
pdf("playground/simulations/mcmc_right_prior_prediction.pdf")
boxplot(ans, main = "Distribution Relative\nPrediction Scores")
dev.off()
