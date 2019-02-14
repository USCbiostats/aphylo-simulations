rm(list = ls())
load("simulations/02-gold-standard/data_and_functions.rda")

library(aphylo)
library(dplyr)
library(magrittr)
library(xtable)

# Setting xtable options
options(xtable.booktabs = TRUE)
options(xtable.include.rownames = FALSE)
options(xtable.sanitize.text.function = function(x) x)

tab <- lapply(dat, "[[", "tip.annotation")
tab <- data.frame(
  tree = rep.int(1:length(dat), sapply(tab, length)),
  anno = unname(do.call(rbind, tab))
)

# CLassifying by size
tab %<>% group_by(tree) %>%
  mutate(
    size = n()
  )
  
quant <- quantile(tab$size, c(.25, .5, .75))
tab %<>% mutate(
  size_bracket = if_else(
    size < quant[1], 1L, 
    if_else(
      size < quant[2], 2L,
      if_else(
        size < quant[3], 3L,
        4L
      )
    ))
)

# Proportions
ans <- tab %>% group_by(tree) %>%
  mutate(
    prop_ones  = sum(anno)/n()
  ) %>%
  ungroup %>%
  group_by(size_bracket) %>%
summarize(
  p025 = quantile(prop_ones, .025),
  p500  = quantile(prop_ones, .5),
  p975 = quantile(prop_ones, .975)
) %>% mutate(
  size_bracket = recode(
    size_bracket,
    "1" = sprintf("(   0,% 4i)", quant[1]),
    "2" = sprintf("[% 4i,% 4i)", quant[1], quant[2]),
    "3" = sprintf("[% 4i,% 4i)", quant[2], quant[3]),
    "4" = sprintf("[% 4i, Inf)", quant[3])
  ) 
) %>%
  mutate(
    size_bracket = paste0("{", size_bracket, "}")
  ) %>%
  rename(
    "Size of the tree" = size_bracket,
    "0.025" = p025,
    "0.50" = p500,
    "0.975" = p975
    ) %T>% print

# Storing the data
xtable(ans) %>%
  `caption<-`("Distribution of ones by tree size") %>%
  `label<-`("fig:tree-size") %>%
  print %>%
  cat(file="tables/tree-size-annotations.tex", sep="\n")
