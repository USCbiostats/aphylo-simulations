
library(dplyr)
library(tidyr)
library(magrittr)
library(aphylo)

source("global-paths.r")

ntrees <- 10

# Reading the list of candidate functions and keeping the trees only
candidate_functions <- readr::read_csv("data/candidate_functions.csv")
candidate_trees     <- unique(candidate_functions$substring)

# Reading panther trees
trees <- lapply(
  sprintf("%s/%s/tree.tree", PANTHER_PATH, candidate_trees[1:ntrees]),
  read_panther
  )

# Preserving the UniProtKB id only
for (i in seq_along(trees))
  trees[[i]]$tree$tip.label <- gsub(".+UniProtKB", "UniProtKB",
                                    trees[[i]]$tree$tip.label)

# Reading the true annotations. We wil merge these with the 
# trees that we will be using.
annotations <- readr::read_delim(
  "data-raw/true_annotations", delim = ";",
  col_names = TRUE)

annotations <- annotations %>% 
  right_join(candidate_functions) %>%
  select(substring, primary_ext_acc, term, qualifier) %>%
  mutate(primary_ext_acc = gsub(".+UniProtKB", "UniProtKB", primary_ext_acc)) 

# Creating aphylo objets
atrees <- vector("list", length(trees))
for (i in seq_len(ntrees)) {
  
  # Gathering the corresponding annotations
  a <- filter(annotations, substring == candidate_trees[i]) %>%
    mutate(
      state = case_when(
        is.na(qualifier) ~ 1L,
        qualifier == "NOT" ~ 0L,
        TRUE ~ 1L
      )
    )
  
  # Reshaping wide
  a <- a %>% 
    select(-substring, -qualifier) %>%
    # Fixing this WEIRD
    # A tibble: 2 x 5
    #   substring primary_ext_acc  term       qualifier      state
    #   <chr>     <chr>            <chr>      <chr>          <int>
    # 1 PTHR10788 UniProtKB=Q0WUI9 GO:0003825 CONTRIBUTES_TO     1
    # 2 PTHR10788 UniProtKB=Q0WUI9 GO:0003825 NOT                0
    group_by(primary_ext_acc, term) %>%
    summarize(state = min(state)) %>%
    spread(term, state)
  
  # Sorting according to the ith tree
  ord <- trees[[i]]$tree$tip.label
  ord <- tibble(id = ord)
  ord <- ord %>%
    left_join(a, by = c("id" = "primary_ext_acc")) %>%
    select(-id)
  
  ord <- as.matrix(ord)
  ord[is.na(ord)] <- 9L
  
  atrees[[i]] <- new_aphylo(
    tip.annotation = ord,
    tree = trees[[i]]$tree
    )
  
}

