library(aphylo)
library(data.table)

# Loading the model parameters
model <- readRDS("novel-predictions/parameter_estimates_mcmc_partially_annotated_no_prior.rds")

# Loading the ADENOSINE DEAMINASE (PTHR11409)
tree  <- read_panther("data-raw/PTHR11409.tree")
tree$tree$tip.label <- gsub(".+UniProtKB", "UniProtKB", tree$tree$tip.label)

# Loading annotations
funs  <- data.table::fread("data-raw/true_annotations")
funs  <- funs[substring == "PTHR11409"]
funs[, label := gsub(".+UniPro", "UniPro", primary_ext_acc)]
funs[is.na(qualifier)  , state := 1L]
funs[qualifier == "NOT", state := 0L] 
funs[is.na(state)      , state := 1L]

funs <- dcast(funs, label ~ term, value.var = "state")
funs <- as.data.frame(funs, stringsAsFactors = FALSE)

# Merging with dataset
tipnames <- data.frame(label = tree$tree$tip.label, stringsAsFactors = FALSE)
funs     <- merge(tipnames, funs, by="label", all.x = TRUE)

for (i. in 1:ncol(funs))
  funs[is.na(funs[[i.]]), i.] <- 9L

# Matching the sorting
funs <- funs[match(tree$tree$tip.label, funs$label), ]
all(funs$label == tree$tree$tip.label)

# Doing the same with the type of node
ntype <- data.frame(
  lab  = rownames(tree$internal_nodes_annotations),
  dupl = as.integer(!tree$internal_nodes_annotations$duplication),
  stringsAsFactors = FALSE
)
ntype <- ntype[match(tree$tree$node.label, rownames(tree$internal_nodes_annotations)),]

stopifnot(all(tree$tree$node.label == ntype$lab))


tree <- new_aphylo(
  tree$tree, tip.annotation = funs[,-1],
  node.type = ntype$dupl
  )

summary(tree)

tree <- tree[c("GO:0004000", "GO:0006154", "GO:0005615")]

samp <- do.call(rbind, model$hist)
idx  <- sample.int(nrow(samp), 1000, TRUE)

pr2 <- predict_pre_order(model, newdata = tree, params = samp[idx[1],])/1000
for (i in 2:length(idx)) {
  pr2 <- pr2 + predict_pre_order(model, newdata = tree, params = samp[idx[i],])/1000
}

pr <- predict(model, newdata = tree)

colnames(pr2) <- paste("Pred. MCMC", colnames(pr))
colnames(pr) <- paste("Pred.", colnames(pr))

tree$tip.annotation <- cbind(tree$tip.annotation, pr[1:Ntip(tree),])
tree$node.annotation <- cbind(tree$node.annotation, pr[-c(1:Ntip(tree)),])
tree$tip.annotation <- cbind(tree$tip.annotation, pr2[1:Ntip(tree),])
tree$node.annotation <- cbind(tree$node.annotation, pr2[-c(1:Ntip(tree)),])

saveRDS(tree, "novel-predictions/PTHR11409.rds")


