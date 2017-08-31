# rm(list = ls())

library(aphylo)

load("simulations/data_and_functions.rda")
load("simulations/mcmc_right_prior_estimates_sample.rda")

N <- 1000

library(parallel)

# Actual proportion of annotations
obs_anns <- lapply(lapply(dat[1:N], "[[", "annotations"), table)
obs_anns <- do.call(rbind, lapply(obs_anns, function(x) {
  cbind(`0` = x["0"], `1` = x["1"])
}))
obs_anns <- as.matrix(obs_anns[complete.cases(obs_anns),])
obs_anns <- colSums(obs_anns)/sum(obs_anns)
prop_of_1s <- obs_anns[2]

cl <- makeForkCluster(10)

# Prediction scores
pred_scores <- parLapply(cl, 1:N , function(i) { 
  try(prediction_score(mcmc_right_prior_estimates_sample[[i]], expected = dat[[i]]$annotations,
                       alpha = prop_of_1s))
})

stopCluster(cl)
save(pred_scores, file = "simulations/mcmc_right_prior_prediction.rda")

ans <- lapply(pred_scores, function(x) {
  tryCatch(cbind(obs = x[["obs"]], rand = x[["random"]]) / x[["worse"]], error = function(e) NULL)
  })
ans <- do.call(rbind, ans)

dimnames(ans) <- list(NULL, c("Model Predictions", "Random Predictions"))


graphics.off()
pdf("simulations/mcmc_right_prior_prediction.pdf")
png("simulations/mcmc_right_prior_prediction.png")
boxplot(ans, main = "Distribution Relative\nPrediction Scores",
        ylab = "Relative Prediction Score\n(0 is perfect prediction)",
        sub  = sprintf(
          "Subsample of 1,000 simulated annotations using PANTHER.\nCompared agains random annotations based on a Bernoulli(%.2f)",
          prop_of_1s)
)
dev.off()
