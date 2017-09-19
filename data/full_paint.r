rm(list = ls())
options(stringsAsFactors = FALSE)

library(aphylo)
library(parallel)
library(readr)
library(stringr)
library(dplyr)


# Reading PAINT annotations ----------------------------------------------------
paint_annotations <- readr::read_tsv(
  "../fullgo_paint081517/Pthr_GO.tsv", n_max=1e7, 
  col_names = c("long_seq_id", "go_id", "qualifier", "evidence_code", "alt_id", "go_ref")
  )

# Subsetting only experimental annotations
experimental_codes <- c("IDA", "EXP", "IPI", "IMP", "IGI", "IEP")
paint_annotations <- subset(paint_annotations, evidence_code %in% experimental_codes)
paint_annotations <- subset(paint_annotations, select = c(long_seq_id, go_id, qualifier))

# Only keep annotations that are yes/no type
paint_annotations <- subset(paint_annotations, !is.na(qualifier))

# Relabeling the annotations
paint_annotations$binary_qualifier <- with(
  paint_annotations, 
  ifelse(grepl("^[^N]", qualifier), 1L, 0L)
)


# Reading trees ----------------------------------------------------------------
tree_fn <- list.files("../PANTHER11.1/books", pattern = "PTHR[0-9]",full.names = TRUE)
tree_fn <- paste0(tree_fn,"/tree.tree")

cl <- makeForkCluster(10)

trees <- parLapply(cl, tree_fn, read_panther)

# Merging the data -------------------------------------------------------------

# First, we need to remove te A[0-9]+ values from the tree ids
leafs <- lapply(trees, function(x) x$tree[["tip.label"]])
leafs <- parLapply(cl, leafs, stringr::str_replace, pattern = "^[a-zA-Z0-9]+[:]",
                replacement = "")
leafs <- parLapply(cl, leafs, function(l)
  data.frame(long_seq_id = l, stringsAsFactors = FALSE)
  )

# Merging annotations
annotations <- parLapply(cl, leafs, function(l) {
  
  # Merging the data for the l-th tree
  unique(dplyr::left_join(l, paint_annotations, by = "long_seq_id"))
})


annotations_tab <- parLapply(cl, seq_along(annotations), function(i, annotations) {
  data.frame(
    tree_id = i,
    n_leafs = length(unique(annotations[[i]]$long_seq_id)),
    annotations[[i]][,c("go_id", "binary_qualifier")]
    , stringsAsFactors = FALSE)
}, annotations=annotations)

annotations_tab <- annotations_tab %>% bind_rows %>%
  filter(!is.na(go_id))

annotations_tab <- 
  annotations_tab %>% group_by(tree_id, go_id) %>%
  summarize(
    no             = sum(binary_qualifier == 0),
    yes            = sum(binary_qualifier == 1),
    prop_annotated = n()/n_leafs[1],
    nleafs         = n_leafs[1]
  ) %>% filter(no > 0 & yes > 0)

# Possible candidates
# tree, go term, no,yes/total leafs
# 453 GO:0003825 5,3/211
# 453 GO:0004805 6,2/211

stopCluster(cl)
