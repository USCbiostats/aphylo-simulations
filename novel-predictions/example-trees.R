library(aphylo)
library(data.table)

files <- list.files(pattern = "^(mcmc|mle).+\\.rds", full.names = FALSE)
estimates <- lapply(files, readRDS)

# Getting the tree names -------------------------------------------------------
tree_names <- sapply(estimates[[4]]$dat, function(i) {
  ids <- which(i$tip.annotatio[, 1] %in% c(0,1))[1]
  i$tree$tip.label[ids]
  })
tree_names <- data.table(UniProtKB = tree_names, ord = 1:length(tree_names))

annotations         <- data.table::fread("../data-raw/true_annotations")
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

ps <- lapply(estimates, prediction_score, loo=TRUE)
aucs <- sapply(ps, function(i) {
  mean(sapply(i, function(j) j$auc$auc))
})
pss  <- sapply(ps, function(i) {
  tmp <- t(sapply(i, function(j) with(j, c(obs, random))))
  colMeans(tmp, na.rm = TRUE)
})

# View(cbind(sapply(ps[[4]], function(i) i$auc$auc), Ntip(estimates[[4]])))

op <- par(xpd = NA)

# Example 1 --------------------------------------------------------------------
graphics.off()
pdf("example-trees-good1-loo.pdf", width = 6, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 11,
  nsamples   = 500,
  prop       = .2,
  cex        = .3,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = TRUE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[11])
  )
dev.off()

graphics.off()
pdf("example-trees-good1.pdf", width = 6, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 11,
  nsamples   = 500,
  prop       = .2,
  cex        = .3,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = TRUE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[11])
)
dev.off()

# Example 2 --------------------------------------------------------------------
graphics.off()
pdf("example-trees-good2-loo.pdf", width = 6, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 68,
  nsamples   = 500,
  prop       = .2,
  cex        = .3,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = TRUE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[68])
)
dev.off()

graphics.off()
pdf("example-trees-good2.pdf", width = 6, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 68,
  nsamples   = 500,
  prop       = .2,
  cex        = .3,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = FALSE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[68])
)
dev.off()

# Example 3 --------------------------------------------------------------------
graphics.off()
pdf("example-trees-bad1-loo.pdf", width = 6, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 124,
  nsamples   = 500,
  prop       = .2,
  cex        = .7,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = TRUE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[24])
)
dev.off()

graphics.off()
pdf("example-trees-bad1.pdf", width = 6, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 124,
  nsamples   = 500,
  prop       = .2,
  cex        = .7,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = FALSE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[24])
)
dev.off()

# Example 4 --------------------------------------------------------------------
graphics.off()
pdf("example-trees-bad2-loo.pdf", width = 6, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 110,
  nsamples   = 500,
  prop       = .2,
  cex        = .7,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = TRUE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[10])
)
dev.off()

graphics.off()
pdf("example-trees-bad2.pdf", width = 6, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 110,
  nsamples   = 500,
  prop       = .2,
  cex        = .7,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = FALSE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[10])
)
dev.off()


par(op)