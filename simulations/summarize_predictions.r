summarize_predictions <- function(x., dat.) {
  
  N <- length(x.)
  
  # Actual proportion of annotations
  obs_anns <- lapply(lapply(dat., "[[", "tip.annotation"), table)
  obs_anns <- do.call(rbind, lapply(obs_anns, function(x) {
    cbind(`0` = x["0"], `1` = x["1"])
  }))
  obs_anns <- as.matrix(obs_anns[complete.cases(obs_anns),])
  obs_anns <- colSums(obs_anns)/sum(obs_anns)
  prop_of_1s <- obs_anns[2]
  
  message("Adjusting etas...")
  x. <- lapply(1:N, function(i) {
    x.[[i]]$par <- c(x.[[i]]$par, eta0 = .5, eta1=.5)
    x.[[i]]
  })
  
  # Prediction scores
  message("Computing prediction scores...")
  structure(
    parallel::mclapply(1:N , function(i) { 
      try(prediction_score(
        x.[[i]],
        expected = with(dat.[[i]], rbind(tip.annotation, node.annotation))))
    }),
    prop_of_1s = prop_of_1s
  )
}