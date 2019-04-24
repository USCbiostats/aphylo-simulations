library(aphylo)
library(dplyr)
library(magrittr)

# Loading functions
source("simulations/00-bias-functions.r")

# Reading the data and calculating AUCs
cross_validation <- readRDS("cv/cross_validation.rds")

cross_validation <- lapply(cross_validation, function(i) {
  i$counts <- aphylo:::fast_table_using_labels(i$expected, c(0, 1, 9)) %>%
    t %>% set_colnames(c(0, 1, 9))
  
  i$expected[i$expected == 9] <- NA
  i
})

aucs <- lapply(cross_validation, "[[", "auc") %>%
  lapply("[", c("auc", "n_used")) %>%
  lapply(unlist, recursive=TRUE) %>%
  do.call(rbind, .) %>%
  as_tibble %>%
  mutate(
    tree = names(cross_validation),
    term = sapply(cross_validation, function(i) colnames(i$expected))
  ) %>%
  arrange(desc(auc))

prescore <- lapply(cross_validation, "[[", "pscore") %>%
  lapply(function(s) {
    with(
      s,
      tibble(
        pscore_obs    = obs[1]/worse,
        pscore_rand   = random/worse,
        ann_expected  = list(expected),
        ann_predicted = list(predicted)
        )
      )
  }) %>%
  bind_rows() %>%
  mutate(tree=names(cross_validation))

# Getting the phylogenetic trees -----------------------------------------------

trees <- lapply(cross_validation, "[[", "estimates") %>%
  lapply("[[", "dat") %>%
  tibble(
    tree = names(.),
    aphylo = .
  ) %>%
  mutate(
    term = sapply(aphylo, function(a) colnames(a$tip.annotation))
  )

# Putting all the information together -----------------------------------------
cross_validation_tab <- cross_validation %>%
  lapply("[[", "counts") %>%
  do.call(rbind, .) %>%
  as_tibble() %>%
  mutate(
    tree = names(cross_validation),
    annotated = (`0` + `1`) #,
    # annotated = annotated #/(annotated + `9`)
  ) %>%
  right_join(aucs, by="tree") %>%
  right_join(prescore, by="tree") %>%
  right_join(trees, by = "tree") %>%
  mutate(
    prop0s     = `0`/(`1` + `0`),
    balance    = abs(.5 - prop0s),
    prop0s_tag = interval_tags(prop0s, c(0, .25, .5, .75, 1))
  ) %>%
  arrange(desc(n_used), desc(balance))

# Making sure that this is giving the right information
stopifnot(all(cross_validation_tab$term.x == cross_validation_tab$term.y ))

cross_validation_tab <- cross_validation_tab %>%
  select(-term.x) %>%
  rename(term = term.y)

saveRDS(cross_validation_tab, "cv/cross_validation_tab.rds", compress = FALSE)
cross_validation_tab %>%
  select(-aphylo, -ann_expected, -ann_predicted) %>%
  readr::write_csv(path = "cv/cross_validation_tab.csv")
