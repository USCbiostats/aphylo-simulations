library(aphylo)

trees <- readRDS("data/candidate_trees.rds")
trees <- do.call(c, trees)
trees <- unlist(lapply(trees, function(d) {
  lapply(1:Nann(d), function(i) d[,i])
}), recursive = FALSE)
trees <- do.call(c, trees)


# Estimating empirical bayes using map
ALPHAS <- c(psi0 = 2, psi1 = 2, mu_d0 = 18, mu_d1 = 10, mu_s0 = 2, mu_s1 = 2, Pi = 2)*3/4
BETAS  <- c(psi0 = 18, psi1 = 18, mu_d0 = 2, mu_d1 = 10, mu_s0 = 18, mu_s1 = 18, Pi = 18)*3/4

bpriors <- bprior(shape1 = ALPHAS, shape2 = BETAS)


map_estimates <- vector("list", Ntrees(trees))
for (i in seq_along(trees)) {
  map_estimates[[i]] <- aphylo_mle(trees[[i]] ~ psi + mu_d + mu_s + Pi, priors = bpriors)
  message("Tree ",i, " of ", Ntrees(trees), " done.")
}

map_estimates <- do.call(rbind, lapply(map_estimates, coef))
ab <- vector("list", ncol(map_estimates))
for (i in 1:ncol(map_estimates)) {
  ab[[i]] <- MASS::fitdistr(
    map_estimates[, i], densfun = "beta",
    start = list(shape1 = ALPHAS[i], shape2 = BETAS[i]),
    lower = .001
  )
}

op <- par(mfrow = c(3, 3))
for (i in 1:length(ab))  {
  curve(dbeta(x, ab[[i]]$estimate[1], ab[[i]]$estimate[2]),
        main = names(ALPHAS)[i], sub = sprintf(
          "MLE: %.4f, Prior: %.4f",
          with(ab[[i]], estimate[1]/sum(estimate)),
          ALPHAS[i]/(ALPHAS[i] + BETAS[i])
        ),
        xlab = "Probability",
        ylab = "Density"
        )
  
  curve(dbeta(x, ALPHAS[i], BETAS[i]),
        main = names(ALPHAS)[i], add = TRUE, col = "red")
}
plot.new()
plot.window(c(0,1),c(0,1))
legend(
  "center",
  cex = 1,
  col = c("black", "red"),
  legend = c(
    "Density based on MLE\nestimates of alpha and beta",
    "Density based on priors\nof alpha and beta"
    )
)

par(op)