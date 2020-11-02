# library(dplyr)
# library(tidyr)
library(magrittr)
library(data.table)
library(aphylo)
# library(slurmR)
# 
source("global-paths.r")

# Reading trees
trees <- readRDS("data/panther_trees.rds")

# Preserving the UniProtKB id only
for (i in seq_along(trees))
  trees[[i]]$tree$tip.label <- gsub(".+UniProtKB", "UniProtKB",
                                    trees[[i]]$tree$tip.label)

# Reading the true annotations. We wil merge these with the 
# trees that we will be using.
annotations         <- data.table::fread("data-raw/true_annotations.gz")
candidate_functions <- data.table::fread("data/candidate_functions.csv")
candidate_trees     <- unique(candidate_functions$substring)

# LEFT JOINT candidate_functions to annotations
annotations <- candidate_functions[annotations, on=c("substring", "term")]
annotations <- annotations[, list(substring, primary_ext_acc, term, qualifier)]
annotations[, primary_ext_acc := gsub(".+UniProtKB", "UniProtKB", primary_ext_acc)]

# Checking competing annotations -----------------------------------------------
annotations[, state := fifelse(is.na(qualifier), 1L, fifelse(qualifier == "NOT", 0L, 1L))]

annotations[, table(qualifier, state)]
# qualifier2
# qualifier             0      1
#                       0 393285
# COLOCALIZES_WITH      0   2694
# CONTRIBUTES_TO        0   2379
# NOT                3092      0

annotations[, off := diff(range(state, na.rm = TRUE)), by = .(primary_ext_acc, term)]
annotations[, table(off)]

# Dropping negati
# off
#      0     1 
# 401424    26
annotations <- annotations[off == 0]
annotations[, c("off") := NULL]

# [2020-11-01] George: PTHR11706 has wrong negative annotations on GO:0010039
#  so we will remove that family+annotation from the analysis.
annotations <- annotations[
  (substring != "PTHR11706" & term != "GO:0010039") |
    (substring == "PTHR11706" & term != "GO:0010039")  
  ]


# Creating aphylo objets
atrees <- vector("list", length(trees))
for (i in seq_along(atrees)) {
  
  # Gathering the corresponding annotations
  a <- annotations[substring == candidate_trees[i]]
  a <- a[, .(state, primary_ext_acc, term)]
  a <- dcast(a, primary_ext_acc ~ term, value.var = "state")

  # Sorting according to the ith tree
  ord <- trees[[i]]$tree$tip.label
  ord <- data.table(id = ord)
  ord <- a[ord, on = c("primary_ext_acc" = "id")]
  ord[, primary_ext_acc := NULL]
  
  # NAs into 9s
  ord <- as.matrix(ord)
  ord[is.na(ord)] <- 9L

  # Now, we match the annotations!
  labs  <- with(trees[[i]]$tree, c(node.label))
  labs  <- data.table(id = labs)
  
  types    <- trees[[i]]$internal_nodes_annotations
  types$id <- rownames(types)
  types    <- data.table(types)
  
  types <- types[labs, on = "id"]
  types <- types[, as.integer(!duplication)]
  types[is.na(types)] <- 0L
  
    
  atrees[[i]] <- new_aphylo(
    tip.annotation = ord,
    tree           = trees[[i]]$tree,
    node.type      = types
    )
  
  if (!(i %% 20))
    message("Tree ", i, " done...")
}

names(atrees) <- candidate_trees
saveRDS(atrees, file = "data/candidate_trees.rds")

# Proportion of annotations, zeros, and ones -----------------------------------
atrees <- readRDS("data/candidate_trees.rds")

candidate_trees <- lapply(names(atrees), function(t.) {
  
  anns <- atrees[[t.]]$tip.annotation
  anns <- apply(anns + 1, 2, tabulate, nbins = 10)[c(1,2,10),]
  
  data.table(
    tree = t.,
    name = colnames(anns),
    zero = anns[1, ],
    one  = anns[2, ],
    miss = anns[3, ]
  )
  
})

candidate_trees <- rbindlist(candidate_trees)
candidate_trees <- candidate_trees[one > 1 & zero > 1]

atrees <- atrees[unique(candidate_trees$tree)]
filtered_trees <- structure(
  vector("list", nrow(candidate_trees)),
  names = candidate_trees$tree
)
for (i in seq_len(nrow(candidate_trees))) {
  
  t. <- candidate_trees$tree[i]
  f. <- candidate_trees$name[i]
  
  filtered_trees[[i]] <- atrees[[t.]][,f.]
  
}

saveRDS(filtered_trees, file = "data/candidate_trees.rds")



