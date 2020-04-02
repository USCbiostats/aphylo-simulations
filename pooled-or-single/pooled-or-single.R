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
ps_pooled_unif <- sapply(ps[[1]], "[[", "obs")
ps_pooled_beta <- sapply(ps[[2]], "[[", "obs")
ps_single_unif <- sapply(estimates_disjoint, function(i) i$ps_unif$obs)
ps_single_beta <- sapply(estimates_disjoint, function(i) i$ps_beta$obs)

ps_obs <- cbind(
  `(a) Pooled-data\n(no prior)` = ps_pooled_unif,
  `(b) Pooled-data\n(beta prior)` = ps_pooled_beta,
  `(c) One-at-a-time\n(no prior)` = ps_single_unif,
  `(d) One-at-a-time\n(beta prior)` = ps_single_beta
)

legends <- rbind(
  t.test(ps_obs[,1], ps_obs[,2], paired = TRUE)$conf.int,
  t.test(ps_obs[,1], ps_obs[,3], paired = TRUE)$conf.int,
  t.test(ps_obs[,1], ps_obs[,4], paired = TRUE)$conf.int
)

legends <- sprintf(
  "(%s) - (%s): [% 4.2f, % 4.2f] %s",
  "a", letters[2:4],
  legends[,1], legends[,2],
  ifelse(sign(legends[,1]) * sign(legends[,2]) > 0, "*", "") 
)
legends <- c(
  "95% C.I. paired t-test",
  legends,
  "* indicates significant difference."
  )

graphics.off()
pdf("pooled-or-single/pooled-or-single.pdf", width = 7, height = 5)
op <- par(mar = par()$mar * c(1, 2, .5,1))
boxplot(
  ps_obs[,4:1], horizontal = TRUE, las = 1,
  xlab = "Prediction score (Leave-one-out)",
  border = gray(.4),
  col    = "gray"
  )

legend(
  "bottomleft", legend = legends,
  bty = "n",
  cex = .75
  )

par(op)
dev.off()
