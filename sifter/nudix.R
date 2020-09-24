library(aphylo)
library(data.table)

model_aphylo <- readRDS("novel-predictions/mcmc_partially_annotated_no_prior.rds")
model_aphylo <- window(model_aphylo, start=5000)

source("sifter/read_sifter.R")

# read.tree("sifter/nudix/PF00293.tree")

tree <- read_nhx("sifter/nudix/PF00293.tree")

# Finding the corresponding annotations
ann  <- fread("sifter/nudix/nudix-ann.tsv")
go <- strsplit(ann$goterms, split=",")
ann <- ann[rep(1:.N, sapply(go, length))]
ann$terms <- unlist(go)
ann[, goterms := NULL]

# Looking for duplicated labs
labs <- data.table(lab = tree$tree$tip.label, name = gsub("[/].+", "", tree$tree$tip.label))
labs[, n := .N, by = name]
setorder(labs, n)

labs[, table(n)] # 10/3693 are duplicated (so we have 5 records)
labs <- merge(
  x = labs,
  y = goa[, .(synm, GO_ID, Aspect, Qualifier)],
  by.x = "name", by.y = "synm")

# How many have experimental annotations?
unique(labs[, .(Aspect, name)])[, table(Aspect)]

dupl <- imputate_duplications(
  tree$tree,
  gsub(".+[_]|[/].+", "", tree$tree$tip.label)
  )

ann[, value:=1L]
ann[, c("moc", "number") := NULL]
anno_i <- dcast(unique(ann), name ~ go, value.var = "value")

# Reading the aphylo tree
ans <- aphylo_from_data_frame(
  tree        = tree$tree,
  annotations = as.data.frame(anno_i),
  types       = data.frame(
    tree$tree$node.label,
    !dupl[-c(1:ape::Ntip(tree$tree))]
  )
)

# Predictions -----------------------------------------------------------------

ids <- which(rowSums(ans$tip.annotation) != Nann(ans) * 9)
pred <- predict(model_aphylo, newdata = ans, ids = ids, loo = TRUE)

pred <- cbind(
  predicted = as.vector(pred[ids,]),
  labels    = as.vector(ans$tip.annotation[ids,])
)

# This is the weird bit
pred[pred==9] <- 0

# Computing AUCs
pred_auc <- auc(pred = pred[,"predicted"], labels = pred[,"labels"])
pred_auc
# Number of observations     : 12656
# Area Under The Curve (AUC) : 0.70
# Rates can be accessed via the $ operator.
plot(pred_auc)
