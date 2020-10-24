library(aphylo)
library(data.table)

estimates <- readRDS("novel-predictions/mcmc_partially_annotated_no_prior.rds")

ids <- lapply(estimates$dat, function(i) {
  which(i$tip.annotation != 9)
})

predlll <- predict(pred, loo = TRUE, ids = ids)
predlll <- lapply(predlll, function(p) {
  p <- p[!is.na(p[,1]),,drop=FALSE]
  data.table(
    name = rownames(p),
    term = colnames(p),
    score = p[,1]
  )
})

predlll <- rbindlist(predlll)

