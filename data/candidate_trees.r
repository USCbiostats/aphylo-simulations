library(dplyr)
library(tidyr)
library(magrittr)
library(aphylo)
library(sluRm)

source("global-paths.r")

ntrees <- 10

# Reading the list of candidate functions and keeping the trees only
candidate_functions <- readr::read_csv("data/candidate_functions.csv")
candidate_trees     <- unique(candidate_functions$substring)

# Reading panther trees
trees <- Slurm_lapply(
  sprintf("%s/%s/tree.tree", PANTHER_PATH, candidate_trees), 
  read_panther,
  njobs      = 10,
  job_name   = "candidate-trees",
  job_path   = STAGING_PATH,
  sbatch_opt = list(account = "lc_pdt", partition = "thomas")
  )

trees <- Slurm_collect(trees)

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
for (i in seq_along(atrees)) {
  
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
  
  # if (!(i %% 10))
  #   message(sprintf("atree[%i] done", i))
  
}

names(atrees) <- candidate_trees
saveRDS(atrees, file = "data/candidate_trees.rds")

# Proportion of annotations, zeros, and ones -----------------------------------
# candidate_trees <- readRDS("data/candidate_trees.rds")

candidate_trees <- lapply(names(candidate_trees), function(t.) {
  
  anns <- candidate_trees[[t.]]$tip.annotation
  anns <- apply(anns, 2, aphylo:::fast_table_using_labels, ids = c(0,1,9))
  tibble(
    tree = t.,
    name = colnames(anns),
    zero = anns[1, ],
    one  = anns[2, ],
    miss = anns[3, ]
  )
  
}) %>% bind_rows %>%
  filter(one > 0, zero > 0)

# Total annotaions
library(ggplot2)
candidate_trees %>%
  mutate(
    zero_prop = zero/(zero + one + miss),
    one_prop  = one/(zero + one + miss),
    miss_prop  = 1 - (zero_prop + one_prop)
  ) %>%
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


