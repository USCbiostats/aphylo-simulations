library(aphylo)
library(dplyr)
library(magrittr)

# Loading functions
source("simulations/00-bias-functions.r")

# Reading the data and calculating AUCs
cross_validation <- readRDS("cv/cross_validation.rds")

candidate_functions <- readr::read_csv("data/candidate_functions.csv")

cross_validation <- lapply(cross_validation, function(i) {
  i$counts <- aphylo:::fast_table_using_labels(i$estimates$dat$tip.annotation, c(0, 1, 9))  %>%
    t %>% set_colnames(c(0, 1, 9))
  
  i$expected[i$expected == 9] <- NA
  i$term <- colnames(i$expected)
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
  mutate(tree = gsub("[-].+$", "", tree)) %>%
  arrange(desc(auc))

prescore <- lapply(cross_validation, "[[", "pscore") %>%
  lapply(function(s) {
    with(
      s,
      tibble(
        term          = colnames(predicted),
        pscore_model  = obs[1],
        pscore_rand   = random,
        ann_expected  = list(expected),
        ann_predicted = list(predicted)
        )
      )
  }) %>%
  bind_rows() %>%
  mutate(
    tree = names(cross_validation),
    tree = gsub("[-][0-9]+$", "", tree) 
    )

# Getting the phylogenetic trees -----------------------------------------------

trees <- lapply(cross_validation, "[[", "estimates") %>%
  lapply("[[", "dat") %>%
  tibble(
    tree = gsub("[-].+$", "", names(.)),
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
    tree = gsub("[-].+$", "", tree),
    term = sapply(cross_validation, "[[", "term"),
    annotated = (`0` + `1`) #,
    # annotated = annotated #/(annotated + `9`)
  ) %>%
  right_join(aucs, by=c("tree", "term")) %>%
  right_join(prescore, by=c("tree", "term")) %>%
  right_join(trees, by = c("tree", "term")) %>%
  mutate(
    prop0s     = `0`/(`1` + `0`),
    balance    = abs(.5 - prop0s),
    prop0s_tag = interval_tags(prop0s, c(0, .25, .5, .75, 1))
  ) %>%
  arrange(desc(n_used), desc(balance))

saveRDS(cross_validation_tab, "cv/cross_validation_tab.rds", compress = FALSE)
cross_validation_tab %>%
  select(-aphylo, -ann_expected, -ann_predicted) %>%
  readr::write_csv(path = "cv/cross_validation_tab.csv")

# Distribution of Prediction scores --------------------------------------------
library(ggplot2)

set.seed(1231)
cross_validation_tab %>%
  mutate(`2 or more` = if_else(`1` > 1 & `0` > 1, "Yes", "No")) %>%
  ggplot(aes(x = pscore_model, y = pscore_rand, color = balance, size=`2 or more`)) +
  geom_jitter(width = .01, height = .01, alpha = .5) +
  geom_abline() +
  labs(x = "Model", y = "Random", color = "Balance", size = "2 or more") +
  theme_bw(base_family = "serif") +
  scale_color_viridis_c() 
  
ggsave("cv/cross_validation_tab.pdf", width = 7, height = 6)  


# Best and worse predictions ---------------------------------------------------

best_and_worse <- cross_validation_tab %>%
  select(tree, term, `0`, `1`, `9`, pscore_model, pscore_rand, auc, aphylo,
         ann_predicted, ann_expected) %>%
  mutate(
    bad  = pscore_rand < pscore_model,
    diff = abs(pscore_rand - pscore_model)
    ) %>%
  group_by(bad) %>%
  arrange(desc(diff)) %>%
  filter(`0` > 1 | `1` > 1, `0` > 0 & `1` > 0) %>%
  filter(1:n() <= 5) %>%
  rename(
    Tree = tree,
    Term = term,
    `n/a` = `9`,
    `Model PS` = pscore_model,
    `Rand. PS` = pscore_rand,
    AUC = auc
  ) %>%
  mutate_at(vars(`0`, `1`, `n/a`), as.integer) %>%
  ungroup %>%
  select(-bad, -diff)

best_and_worse %>%
  select(-aphylo, -ann_predicted, -ann_expected) %>%
  xtable::xtable(
    caption = paste(
      "Selected cases where the prediction model performed the best and the worse.",
      "It is important to point out that there are no cases in which there were 2",
      "or more annotations of each and the model performed worse than random."),
    label = "tab:top10"
  ) %>%
  xtable::print.xtable(
    file = "cv/cross_validation_top10.tex",
    include.rownames = FALSE,
    booktabs = TRUE
    )


tibble(
  Observed  = best_and_worse$ann_expected[[1]][,1],
  Predicted = best_and_worse$ann_predicted[[1]][,1]
) %>%
  arrange(Predicted) %>%
  xtable::xtable(
    caption = sprintf("Tree: %s, Term: %s", best_and_worse$Tree[1], best_and_worse$Term[1])
  ) %>%
  xtable::print.xtable(
    file = "cv/cross_validation_example_best.tex",
    include.rownames = TRUE
  ) 

tibble(
  Observed  = best_and_worse$ann_expected[[6]][,1],
  Predicted = best_and_worse$ann_predicted[[6]][,1]
) %>%
  arrange(Predicted) %>%
  xtable::xtable(
    caption = sprintf("Tree: %s, Term: %s", best_and_worse$Tree[6], best_and_worse$Term[6])
  ) %>%
  xtable::print.xtable(
    file = "cv/cross_validation_example_worst.tex",
    include.rownames = TRUE
  ) 


graphics.off()
pdf("cv/cross_validation_example_best.pdf", width=6, height=6)
plot(
  best_and_worse$aphylo[[1]], show.node.label = FALSE, show.tip.label = FALSE,
  rect.args = list(border="transparent"),
  main = sprintf("Tree: %s\nTerm: %s", best_and_worse$Tree[1], best_and_worse$Term[1])
  )
dev.off()

pdf("cv/cross_validation_example_worse.pdf", width=6, height=6)
plot(
  best_and_worse$aphylo[[6]], show.node.label = FALSE, show.tip.label = FALSE,
  rect.args = list(border="transparent"),
  main = sprintf("Tree: %s\nTerm: %s", best_and_worse$Tree[6], best_and_worse$Term[6])
)
dev.off()

