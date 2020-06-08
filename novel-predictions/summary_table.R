library(aphylo)

# Reading the results
files <- list.files(
  path       = "novel-predictions/",
  pattern    = "^(mcmc|mle).+partially.+\\.rds",
  full.names = TRUE
  )
fnames <- list.files(
  path       = "novel-predictions/",
  pattern    = "^(mcmc|mle).+partially.+\\.rds",
  full.names = FALSE
)

estimates <- lapply(files, readRDS)

# Getting these as names
fnames <- gsub("(mcmc|mle)", "\\U\\1", fnames, perl = TRUE)
fnames <- gsub("_", " ", fnames)
fnames <- gsub("[.]rds", "", fnames)
names(estimates) <- fnames

# Splitting chains
for (i in which(grepl("^MCMC", fnames))) {
  estimates[[i]] <- window(
    x     = estimates[[i]],
    start = floor(coda::niter(estimates[[i]]$hist) / 2)
  )
}

theta <- do.call(cbind, lapply(estimates, coef))
vars  <- do.call(cbind, lapply(lapply(estimates, vcov), diag))
vars[] <- sqrt(vars)

# Adding information
ps <- parallel::mclapply(
  estimates,
  prediction_score,
  loo = TRUE, mc.cores = 4L
)

tab  <- matrix(sprintf("%.2f", theta), ncol = ncol(theta))
dimnames(tab) <- list(rownames(theta), paste0("(", 1:ncol(tab), ")"))

varnames <- c(
  psi0  = "$\\psi_{01}$",
  psi1  = "$\\psi_{10}$",
  mu_d0 = "$\\mu_{d01}$",
  mu_d1 = "$\\mu_{d10}$",
  mu_s0 = "$\\mu_{s01}$",
  mu_s1 = "$\\mu_{s10}$",
  Pi    = "$\\pi$"
)

rownames(tab) <- varnames[rownames(tab)]

aucs <- sapply(ps, function(ps.) {
  a <- do.call(rbind, lapply(ps., function(p) cbind(p$expected, p$predicted)))
  auc(a[,2],labels = a[,1])$auc
})

aucs10 <- sapply(ps, function(ps.) {
  a <- lapply(ps., function(p) cbind(p$expected, p$predicted))
  a <- lapply(a, function(a.) a.[!is.na(a.[,2]),])
  a <- a[sapply(a, function(i) sum(1-i[,2])) >= 10]
  print(length(a))
  a <- do.call(rbind, a)
  # a <- a[!is.na(a[,2])]
  auc(a[,2],labels = a[,1])$auc
})

maes <- sapply(ps, function(ps.) {
  a <- do.call(rbind, lapply(ps., function(p) cbind(p$expected, p$predicted)))
  a <- a[!is.na(a[,2]),]
  mean(abs(a[,1] - a[,2]))
})

# aucs <- sapply(ps, function(i) {
#   mean(sapply(i, function(j) j$auc$auc))
# })
# pss  <- sapply(ps, function(i) {
#   tmp <- t(sapply(i, function(j) with(j, c(obs, random))))
#   colMeans(tmp, na.rm = TRUE)
# })


tab <- rbind(
  tab,
  `Tree count` = sapply(estimates, Ntrees),
  Method       = ifelse(grepl("mcmc", files), "MCMC", "MLE"),
  Prior        = ifelse(grepl("no_prior", files), "Uniform", "Beta"),
  # Inferred     = ifelse(grepl("fully", files), "Yes", "No"),
  AUC          = sprintf("%.2f", aucs),
  `MAE`        = sprintf("%.2f", maes)
)

is_map <- tab["Method",] == "MLE" & tab["Prior",] == "Beta"
is_mle <- tab["Method",] == "MLE" & tab["Prior",] == "Uniform"
tab["Prior", is_mle]  <- "-"
tab["Method", is_map] <- "MAP"

print(
  xtable::xtable(tab),
  sanitize.text.function = function(i) i,
  booktabs = TRUE
)
