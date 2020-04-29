library(aphylo)
library(data.table)

# Listing the mcmc files
files <- list.files(
  pattern = "^(mcmc|mle).+\\.rds",
  path    = "novel-predictions",
  full.names = TRUE
)
estimates <- lapply(files, readRDS)
estimates_disjoint <- readRDS("novel-predictions/estimates_disjoint.rds")


# Getting the tree names -------------------------------------------------------
tree_names <- sapply(estimates[[4]]$dat, function(i) {
  ids <- which(i$tip.annotatio[, 1] %in% c(0,1))[1]
  i$tree$tip.label[ids]
})
tree_names <- data.table(UniProtKB = tree_names, ord = 1:length(tree_names))

annotations         <- data.table::fread("data-raw/true_annotations")
annotations[, UniProtKB := gsub(".+UniProtKB", "UniProtKB", primary_ext_acc)]
annotations <- subset(annotations, select = c(substring, UniProtKB))
annotations <- unique(annotations)

# Duplicates?
table(table(annotations$UniProtKB)) # Should equal to 1

tree_names <- merge(annotations, tree_names, all.x = FALSE, all.y = TRUE, by = "UniProtKB")
tree_names <- as.data.frame(tree_names)[order(tree_names$ord),]$substring

# Creating a copy, we need the full chains for later
estimates0 <- estimates

# Splitting chains
for (i in which(grepl("^mcmc", sapply(estimates, "[[", "method")))) {
  estimates[[i]] <- window(
    x     = estimates[[i]],
    start = (end(estimates[[i]]$hist) - start(estimates[[i]]$hist))/2
  )
}

# Getting which files are MCMC partially annotated
ids_mcmc_partially <- which(grepl("mcmc.+partially", files))

ps <- lapply(estimates[ids_mcmc_partially], prediction_score, loo=TRUE)

# Extracting prediction scores and comparing -----------------------------------

# Joint model
ps_pooled_unif <- sapply(ps[[1]], function(p) p$obs) # "[[", "obs")
ps_pooled_beta <- sapply(ps[[2]], function(p) p$obs) # "[[", "obs")

ps_single_unif <- sapply(estimates_disjoint, function(i) {
  with(i$ps_unif, {
    
    n  <- length(obs.ids)
    a0 <- max(sum(expected[obs.ids,] == 0) - 1, 0)/(n - 1)
    a1 <- max(sum(expected[obs.ids,] == 1) - 1, 0)/(n - 1)
    
    prediction_score(
      x        = predicted[obs.ids,,drop=FALSE],
      expected = expected[obs.ids,,drop=FALSE],
      alpha0   = a0,
      alpha1   = a1
      )
  })$obs
})


ps_single_beta <- sapply(estimates_disjoint, function(i) {
  with(i$ps_beta,{
    
    n  <- length(obs.ids)
    a0 <- max(sum(expected[obs.ids,] == 0) - 1, 0)/(n - 1)
    a1 <- max(sum(expected[obs.ids,] == 1) - 1, 0)/(n - 1)
    
    prediction_score(
      x        = predicted[obs.ids,,drop=FALSE],
      expected = expected[obs.ids,,drop=FALSE],
      alpha0   = a0,
      alpha1   = a1
      )
    })$obs
  })

i <- 4
n <- length(with(ps[[1]][[i]], expected[obs.ids]))
prandom <- 1 - aphylo:::predict_random(
  1, with(ps[[1]][[i]], expected[obs.ids,,drop=FALSE]),
  G_inv = diag(n),
  alpha0 = with(ps[[1]][[i]], (sum(expected[obs.ids,] == 0) - 1)/(n - 1)),
  alpha1 = with(ps[[1]][[i]], (sum(expected[obs.ids,] == 1) - 1)/(n - 1)),
  R      = 5e4
  )/n
mean(prandom)
1 - aphylo:::prediction_score_rand(
  with(ps[[1]][[i]], expected[obs.ids,,drop=FALSE]),
  W = diag(n),
  alpha0 = with(ps[[1]][[i]], (sum(expected[obs.ids,] == 0) - 1)/(n - 1)),
  alpha1 = with(ps[[1]][[i]], (sum(expected[obs.ids,] == 1) - 1)/(n - 1))
)/n

ps_random      <- sapply(ps[[2]], "[[", "random")

ps_random2     <- sapply(estimates_disjoint, function(i) i$ps_beta$random)

ps_obs <- cbind(
  `(random) Random\n(unfair coin toss)`= ps_random,
  `(a) Pooled-data\n(Uniform prior)`   = ps_pooled_unif,
  `(b) Pooled-data\n(Beta prior)`      = ps_pooled_beta,
  `(c) One-at-a-time\n(Uniform prior)` = ps_single_unif,
  `(d) One-at-a-time\n(Beta prior)`    = ps_single_beta
)

legends <- rbind(
  t.test(ps_obs[,2], ps_obs[,1], paired = TRUE)$conf.int,
  t.test(ps_obs[,2], ps_obs[,3], paired = TRUE)$conf.int,
  t.test(ps_obs[,2], ps_obs[,4], paired = TRUE)$conf.int,
  t.test(ps_obs[,2], ps_obs[,5], paired = TRUE)$conf.int
)

legends <- sprintf(
  "(%s) - (%s): [% 4.2f, % 4.2f] %s",
  "a", c("random", letters[2:4]),
  legends[,1], legends[,2],
  ifelse(sign(legends[,1]) * sign(legends[,2]) > 0, "*", "") 
)
legends <- c(
  "95% C.I. paired t-test",
  legends,
  "* indicates significant difference."
  )

graphics.off()
pdf("pooled-or-single/pooled-or-single.pdf", width = 7, height = 7)
op <- par(mar = par()$mar * c(1, 2, .5,1))
boxplot(
  ps_obs[,5:1], horizontal = TRUE, las = 1,
  xlab = "Prediction score (Leave-one-out)",
  border = gray(.4),
  col    = "gray",
  ylim   = c(0, 1)
  )

legend(
  "topleft", legend = legends,
  bty = "n",
  cex = .75
  )

par(op)
dev.off()

tab <- do.call(rbind, lapply(ps[[1]], function(p) {
  table(with(p, expected[obs.ids,]))
}))

# Checking cases in which AUC is either 0/1
aucs_i <- sapply(ps[[1]], function(p) p$auc$auc)
pval_i <- sapply(ps[[1]], function(p) p$pval)
ps_i   <- sapply(ps[[1]], function(p) p$obs)
rand_i <- sapply(ps[[1]], function(p) p$random)

dat_i <- cbind(aucs_i, pval_i, ps_i, rand_i)
View(dat_i[dat_i[,1] > .9 | dat_i[,1] < .1,])

# Creating a table instead looking at pairwise statistics -----------------------

S <- ps_obs
S <- 1-S # Let's make it MAE again...
colnames(S) <- gsub("^[(][a-zA-Z]+[)]\\s+", "", colnames(S))
colnames(S) <- gsub("\\n", " ", colnames(S), fixed = FALSE)

differences <- matrix("-", nrow = ncol(ps_obs), ncol = ncol(ps_obs),
                      dimnames = list(colnames(S), colnames(S)))
for (i in 1:ncol(S)) {
  for (j in 1:ncol(S)) {
    if (i == j)
      next
    
    # Testing the differences
    if (i == 1)
      tmp <- t.test(
        x           = S[,j],
        mu          = 0.5,
        alternative = "two.sided"
        )
    else if (j == 1)
      tmp <- t.test(
        x           = S[,i],
        mu          = 0.5,
        alternative = "two.sided"
      )
    else
      tmp <- t.test(
        x           = S[,i],
        y           = S[,j],
        paired      = TRUE,
        var.equal   = FALSE,
        alternative = "two.sided"
      )
    
    # Adjusting the confint
    if (i == 1)
      tmp$conf.int <- .5 - tmp$conf.int
    else if (j == 1)
      tmp$conf.int <- tmp$conf.int - .5
    
    message(i, " ", j)
    
    differences[i, j] <- sprintf(
      "%s[% 4.2f,% 4.2f]%s",
      if (prod(tmp$conf.int) > 0) {
        if (tmp$conf.int < 0) "\\cellcolor{blue!25}"
        else "\\cellcolor{red!25}"
      } else
        "",
      tmp$conf.int[1], tmp$conf.int[2],
      with(tmp, {
        ifelse(p.value < .01, "$^{***}$",
               ifelse(p.value < .05, "$^{**}$",
                      ifelse(p.value < .1, "$^{*}$", "")))
      })
    )
    
  }
}

differences

# Latexing the table
differences[] <- gsub("([[]|[,])(\\s)", "\\1\\\\hphantom{-}", differences[])
differences <- cbind(
  gsub("^.+[(](.+)[)]$", "\\\\hspace{1mm}\\1", rownames(differences), perl = TRUE),
  differences
)
differences[,1] <- gsub("Uniform", "Unif.", differences[,1], perl = TRUE)

f <- tempfile()
print(
  xtable::xtable(differences), file = f,
  sanitize.text.function = function(i) i,
  booktabs = TRUE, include.rownames = FALSE
  )
tab <- readLines(f)

# Manual updates
tab[grepl("begin{tabular}", tab, fixed = TRUE)]  <- "\\begin{tabular}{r*{5}{m{.15\\linewidth}}}"
tab[grepl("toprule", tab)] <- "\\toprule & & \\multicolumn{2}{c}{Pooled-data} & \\multicolumn{2}{c}{One-at-a-time} & \\"
cat(tab, sep = "\n")
