library(aphylo)

# Preparing the data
proportions <- c(.1, .2, .4, .5, .7, .8, .9)
n <- 20

A <- lapply(proportions, function(p) {
  n1 <- p * n
  matrix(c(
    rep(0, n - n1),
    rep(1, n1)
  ), ncol = 1)
})

sapply(A, mean)
# alphas <- seq(from=.01, to=.99, length.out = 8)

calc_a0 <- function(A) {
  n <- nrow(A)
  max(sum(A == 0) - 1, 0)/(n - 1)
}

calc_a1 <- function(A) {
  n <- nrow(A)
  max(sum(A == 1) - 1, 0)/(n - 1)
}

# Running the experiments
experiments <- vector("list", length(proportions))
for (i in seq_along(experiments)) {
  
  experiments[[i]] <- lapply(proportions, function(a) {
    aphylo:::predict_random(
      1L, A = A[[i]], G_inv = diag(n),
      alpha0 = 1 - a, #calc_a0(A[[i]]),
      alpha1 = a, #calc_a1(A[[i]]),
      R = 5000
      )
  })
  
  experiments[[i]] <- do.call(cbind, experiments[[i]])
  experiments[[i]] <- 1 - experiments[[i]]/n
  colnames(experiments[[i]]) <- sprintf("%.2f", proportions)
  message("Experiment ", i, " done.")
}

op <- par(mfrow = c(4, 2), mar = c(2.5,2,2.5,1), oma = c(4,2,0,0))
for (i in seq_along(experiments)) {
  
  col. <- (which.min(abs(proportions[i] - proportions)) == 1:length(experiments)) + 1
  
  boxplot(
    x    = experiments[[i]],
    ylim = c(0, 1.25),
    col  = c("white", "gray")[col.]
    )
  
  legend("topleft", bty = "n", legend = proportions[i])
  
  abline(h=.5)
  
  # abline(lty = 2, lwd = 2, v = val_i)
}
par(op)


# Experiment 2: Various sizes, looking at the expected value instead, 

ntrees <- 100
sizes  <- 10:100
props  <- c(.1, .5, .9)

graphics.off()
pdf("fig/predscore-analysis.pdf", width = 4, height = 6)
op <- par(mfrow = c(3, 1), mar = c(1, 4, 1, 1), oma = c(6,3,1,1), cex.axis=.8)
set.seed(1231)
for (p in props) {

  A <- replicate(ntrees, {
    cbind(
      as.integer(runif(sample(sizes, 1)) > p)
    )
  }, simplify = FALSE)
  
  alphas <- c(.05, .1, .2, .4, .5, .7, .8, .9, .95)
  experiments <- vector("list", length(alphas))
  for (i in seq_along(alphas)) {
    
    experiments[[i]] <- sapply(A, function(a) {
      prediction_score(cbind(rep(alphas[i], nrow(a))), a)$obs
    })
    message("Repetition ", i, " complete")
  }
  
  experiments <- do.call(cbind, experiments)
  colnames(experiments) <- alphas
  boxplot(
    experiments,
    ylab = sprintf("P(X=1) = %.2f", p), 
    xlab = "",
    border = gray(.4),
    col    = "gray"
    )
  
}
par(op)
title(xlab = expression(Value~of~alpha), ylab = "Prediction score")
dev.off()