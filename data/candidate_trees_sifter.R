library(data.table)
library(aphylo)

trees <- readRDS("data/candidate_trees.rds")
annotations <- fread("data-raw/true_annotations.gz")

annotations[, primary_ext_acc := gsub(".+(?=UniP)", "", primary_ext_acc, perl = TRUE)]

tree_data <- lapply(seq_along(trees), function(i) {
  d   <- trees[[i]]$tip.annotation
  ids <- which(d[,1] != 9)
  data.table(
    term = colnames(d),
    id   = trees[[i]]$tree$tip.label[ids],
    tree = names(trees)[i],
    ann  = d[ids,]
  )
})

tree_data <- rbindlist(tree_data)

# Pasting annotations ----------------------------------------------------------
extended_data <- merge(
  x = unique(tree_data[,.(id,tree)]),
  y = annotations[qualifier != "NOT"],
  by.x = c("id", "tree"),
  by.y = c("primary_ext_acc", "substring"),
  all.x = TRUE,
  all.y = FALSE
)

# Counting how many functions per tree
funs_per_tree <- unique(extended_data[, .(tree, term)])

# SIFTER is OK with between 2 to 8 functions
sifter_relevant_trees <- funs_per_tree[, .(n = .N), by = tree][n %inrange% c(2,10),] 

extended_data <- extended_data[tree %in% sifter_relevant_trees$tree]

# And we should have at least one annotation per protein, so if there's
# any missing is because the annotation was a 0 (which we drop earlier)
extended_data <- extended_data[!is.na(term)]


# Now, how many proteins per tree at this stage?
unique(extended_data[, .(id, tree)])[,  .(n = .N), by = tree]
#         tree n
# 1: PTHR33146 3
# 2: PTHR10961 4
# 3: PTHR10209 7
# 4: PTHR11783 4
# 5: PTHR23342 4
# 6: PTHR10424 4
# > 3+4+7+4+4+4
# [1] 26

