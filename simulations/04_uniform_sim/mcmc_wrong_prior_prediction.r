
library(aphylo)


load("simulations/04_uniform_sim/data_and_functions.rda")
load("simulations/04_uniform_sim/mcmc_wrong_prior_estimates.rda")
source("simulations/summarize_predictions.r")

N <- length(ans_MCMC_wrong_prior)

# Prediction scores
pred_scores <- summarize_predictions(ans_MCMC_wrong_prior, dat, N)
saveRDS(pred_scores, file = "simulations/04_uniform_sim/mcmc_wrong_prior_prediction.rds", compress = FALSE)


ans <- lapply(pred_scores, function(x) {
  tryCatch(cbind(obs = x[["obs"]], rand = x[["random"]]) / x[["worse"]], error = function(e) NULL)
  })
ans <- do.call(rbind, ans)

dimnames(ans) <- list(NULL, c("Model Predictions", "Random Predictions"))


graphics.off()
pdf("simulations/04_uniform_sim/mcmc_wrong_prior_prediction.pdf")
# png("simulations/mcmc_wrong_prior_prediction.png")
boxplot(ans, #main = "Distribution Relative\nPrediction Scores",
        ylab = "Relative Prediction Score (0 is perfect prediction)",
        sub  = sprintf(
          "Random annotations based on a Bernoulli(%.2f)",
          attr(pred_scores, "prop_of_1s")),
        ylim = c(0,.6),outline = FALSE, notch=FALSE
)

# Examples of predictions
set.seed(2)
S <- sample.int(N, 4)
graphics.off()
pdf("simulations/04_uniform_sim/mcmc_wrong_prior_prediction_plots.pdf")
op <- par(mfrow=c(2,2))
for (s in S)
  plot(pred_scores[[s]], main = sprintf("Predicted vs truth tree #%i", s))
par(op)
dev.off()
