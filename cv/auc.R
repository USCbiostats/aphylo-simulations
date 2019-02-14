library(aphylo)
library(dplyr)
library(magrittr)

# Reading the data and calculating AUCs
cross_validation <- readRDS("cv/cross_validation.rds")

cross_validation <- lapply(cross_validation, function(i) {
  i$counts <- aphylo:::fast_table_using_labels(i$expected, c(0, 1, 9)) %>%
    t %>% set_colnames(c(0, 1, 9))
    
  i$expected[i$expected == 9] <- NA
  i
})

dat <- lapply(cross_validation, function(i) {
  with(i, auc(pred = pred_out, expected))
})

aucs <- lapply(dat, "[", c("auc", "n_used")) %>%
  lapply(unlist, recursive=TRUE) %>%
  do.call(rbind, .) %>%
  as_tibble %>%
  mutate(
    tree = names(cross_validation),
    term = sapply(cross_validation, function(i) colnames(i$expected))
    ) %>%
  arrange(desc(auc))

accuracy <- cross_validation %>%
  lapply("[[", "counts") %>%
  do.call(rbind, .) %>%
  as_tibble() %>%
  mutate(
    tree = names(cross_validation),
    annotated = (`0` + `1`),
    annotated = annotated/(annotated + `9`)
  ) %>%
  right_join(aucs, by="tree") %>%
  arrange(desc(auc))
  

dat[[accuracy$tree[9]]] %>%
  plot
  
accuracy %>%
  mutate(prop0s = `0`/(`1` + `0`)) %>%
  filter(prop0s >= .3, prop0s <= .7, n_used >= 8) %>%
  View
