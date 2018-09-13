rm(list = ls())

library(aphylo)


dat <- readRDS("simulations/dgp.rds")
ans_MCMC_right_prior <- readRDS("simulations/04-full-model/mcmc_right_prior_estimates.rds")

source("simulations/summarize_predictions.r")

N <- length(ans_MCMC_right_prior)

# Prediction scores
pred_scores <- summarize_predictions(ans_MCMC_right_prior, dat, N)
saveRDS(pred_scores, file = "simulations/04-full-model/mcmc_right_prior_prediction.rds", compress = FALSE)


ans <- lapply(pred_scores, function(x) {
  tryCatch(cbind(obs = x[["obs"]], rand = x[["random"]]) / x[["worse"]], error = function(e) NULL)
  })
ans <- do.call(rbind, ans)

dimnames(ans) <- list(NULL, c("Model Predictions", "Random Predictions"))


graphics.off()
pdf("simulations/04-full-model/mcmc_right_prior_prediction.pdf")
# png("simulations/mcmc_right_prior_prediction.png")
boxplot(ans, #main = "Distribution Relative\nPrediction Scores",
        ylab = "Relative Prediction Score (0 is perfect prediction)",
        sub  = sprintf(
          "Random annotations based on a Bernoulli(%.2f)",
          attr(pred_scores, "prop_of_1s")),
        ylim = c(0,.6),outline = FALSE, notch=FALSE
)
dev.off()
