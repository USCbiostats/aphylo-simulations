
library(AUC)

aphylo_auci <- function(x) {#, thr = seq(.0001, .9999, length.out = 100L)) {
  
  # It only works for aphylo objects
  if (!inherits(x, "aphylo_prediction_score"))
    return(NA)
  
  dat <-  with(x, cbind(predicted, expected))
  dat[dat[] == 9L] <- NA
  dat <- dat[complete.cases(dat),,drop=FALSE]
  
  if (!length(sum(dat[,2] == 0L)*sum(dat[,2] == 1L)))
    return(NA)
  
  # Calculating AUC alar
  return(try(AUC::auc(AUC::roc(dat[,1], as.factor(dat[,2]))), silent = TRUE))
  
  # # Indices
  # zeros <- which(dat[,2] == 0L)
  # ones  <- which(dat[,2] == 1L)
  # 
  # nzeros <- length(zeros)
  # nones  <- length(ones)
  # 
  # # Computing dots
  # ans <- matrix(ncol=2, nrow=length(thr), dimnames = list(NULL, c("TPR", "FNR")))
  # for (i in seq_along(thr)) {
  #   
  #   pred     <- as.integer(dat[,1] > thr[i])
  #   ans[i, ] <- cbind(
  #     TPR = sum(pred[ones]),
  #     FNR = sum(pred[zeros] != dat[zeros,2])
  #   )
  #   
  # }
  # 
  # ans[,1] <- ans[,1]/nones
  # ans[,2] <- ans[,2]/nzeros
  # 
  # 
  # ans
  
  
}

aphylo_auc <- function(L) {
  
  ans <- parallel::mclapply(L, aphylo_auci, mc.cores = 10L)
  ans[!sapply(ans, inherits, what="numeric")] <- NA
  unlist(ans)
  
}
