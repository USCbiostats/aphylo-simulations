library(aphylo)
library(data.table)

files <- list.files(
  pattern = "^(mcmc|mle).+\\.rds",
  path    = "novel-predictions",
  full.names = TRUE
  )
estimates <- lapply(files, readRDS)

# Getting the tree names -------------------------------------------------------
tree_names <- sapply(estimates[[4]]$dat, function(i) {
  ids <- which(i$tip.annotatio[, 1] %in% c(0,1))[1]
  i$tree$tip.label[ids]
  })
tree_names <- data.table(UniProtKB = tree_names, ord = 1:length(tree_names))

annotations         <- data.table::fread("data-raw/true_annotations")
annotations[, UniProtKB := gsub(".+UniProtKB", "UniProtKB", primary_ext_acc)]

# Temp pretty table
pretty_ann <- annotations[, list(substring, primary_ext_acc, term, qualifier)]
colnames(pretty_ann) <- c("Family", "Id", "GO term", "Qualifier")
pretty_ann <- xtable::xtable(pretty_ann[1:10])
xtable::print.xtable(pretty_ann, booktabs = TRUE)


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
aucs <- sapply(ps, function(i) {
  mean(sapply(i, function(j) j$auc$auc))
})
pss_unif  <- sapply(ps[[1]], "[[", "obs")
pss_beta  <- sapply(ps[[2]], "[[", "obs")

# Is the difference significant?
t.test(pss_unif, pss_beta, paired = TRUE)

View(cbind(
  obs  = sapply(ps[[1]], function(i) i$obs),
  size = Ntip(estimates[[3]]),
  nann = sapply(ps[[1]], function(i) length(i$obs.ids))
))


op0 <- par(mai = rep(0, 4))
graphics.off()
pdf("featured-trees/example-trees-good1.pdf", width = 9, height = 6)
plot(
  x          = estimates[[4]],
  which.tree = 11,
  prop       = .2,
  cex        = .3,
  node.type.size = c(1.25,.75)/20,
  ncores     = 4L,
  loo        = FALSE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[11])
)
dev.off()
par(op0)

op <- par(xpd = NA)

# Example 1 --------------------------------------------------------------------
graphics.off()
svg("featured-trees/example-trees-good1-loo.svg", width = 6, height = 6)
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
svg("featured-trees/example-trees-bad1-loo.svg", width = 6, height = 6)
plot(
  x          = estimates[[3]],
  which.tree = 51, #124,
  nsamples   = 10,
  prop       = .2,
  cex        = .7,
  node.type.size = c(1.25,.75)/30,
  ncores     = 4L,
  loo        = TRUE,
  show.tip.label = FALSE,
  main       = sprintf("%s", tree_names[51])
)
dev.off()


par(op)