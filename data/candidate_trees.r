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
annotations         <- data.table::fread("data-raw/true_annotations")
candidate_functions <- data.table::fread("data/candidate_functions.csv")
candidate_trees     <- unique(candidate_functions$substring)

# LEFT JOINT candidate_functions to annotations
annotations <- candidate_functions[annotations, on=c("substring", "term")]
annotations <- annotations[, list(substring, primary_ext_acc, term, qualifier)]
annotations[, primary_ext_acc := gsub(".+UniProtKB", "UniProtKB", primary_ext_acc)]

# Creating aphylo objets
atrees <- vector("list", length(trees))
for (i in seq_along(atrees)) {
  
  # Gathering the corresponding annotations
  a <- annotations[substring == candidate_trees[i]]
  a[is.na(qualifier)  , state := 1L]
  a[qualifier == "NOT", state := 0L] 
  a[is.na(state)      , state := 1L]
  
  a <- a[, list(state = min(state)), by = c("primary_ext_acc", "term")]
  
  a <- dcast(a, primary_ext_acc ~ term, value.var = "state")
  # Reshaping wide
  # a <-
    # a %>% 
    # select(-substring, -qualifier) %>%
    # # Fixing this WEIRD
    # # A tibble: 2 x 5
    # #   substring primary_ext_acc  term       qualifier      state
    # #   <chr>     <chr>            <chr>      <chr>          <int>
    # # 1 PTHR10788 UniProtKB=Q0WUI9 GO:0003825 CONTRIBUTES_TO     1
    # # 2 PTHR10788 UniProtKB=Q0WUI9 GO:0003825 NOT                0
    # group_by(primary_ext_acc, term) %>%
    # summarize(state = min(state)) %>%
    # spread(term, state)
  
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
saveRDS(atrees, file = "data/candidate_trees2.rds")

# Proportion of annotations, zeros, and ones -----------------------------------
atrees <- readRDS("data/candidate_trees2.rds")

candidate_trees <- lapply(names(atrees), function(t.) {
  
  anns <- atrees[[t.]]$tip.annotation
  anns <- apply(anns, 2, aphylo:::fast_table_using_labels, ids = c(0,1,9))
  
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

saveRDS(filtered_trees, file = "data/candidate_trees2.rds")


# Total annotaions
library(ggplot2)

candidate_trees[, zero_prop := zero/(zero + one + miss)]
candidate_trees[, one_prop  := one/(zero + one + miss)]
candidate_trees[, miss_prop := 1 - (zero_prop + one_prop)]


candidate_trees %>%
  ggplot(aes(x = one_prop, y = zero_prop)) +
  # theme_minimal() +
  # theme(panel.background = element_rect("gray")) +
  geom_hex(bins=20) +
  # scale_fill_distiller() +
  theme_bw() +
  scale_fill_viridis_c() +
  scale_x_log10() +
  scale_y_log10() + 
  xlab("% annotated with 1 (log-scale)") +
  ylab("% annotated with 0 (log-scale)") +
  labs(fill = "Freq") +
  geom_abline(slope = 1, intercept = 0, color="white", lwd=1) 

ggsave(filename = "data/candidate_trees.pdf", width = 7,
       height = 6)


