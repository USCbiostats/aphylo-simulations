

library(ape)
library(aphylo)
library(stringi)

# Functions to read sifter trees
source("sifter/read_sifter.R")

# Reading trees
fn <- list.files("sifter/hundred/", pattern = "nhx", full.names = TRUE)
tn <- gsub("(.+/)([a-zA-Z0-9]+)\\.nhx", "\\2", fn, perl = TRUE)

trees <- parallel::mclapply(fn, read_nhx, mc.cores = 4L)

# How many diplication types
ndupl <- lapply(trees, function(tree) {
  sapply(tree$nhx, function(n) {
    if ("D" %in% names(n) && n["D"] != "N")
      TRUE
    else
      FALSE
  })
})

sapply(ndupl, sum) # No tree here has duplication events

info  <- parallel::mclapply(fn, function(i) {
  x <- readLines(i)
  stringr::str_extract_all(x, "\\[:\\][0-9\\.]+\\[&&NHX:[A-Z][:=][A-Z]\\]")
}, mc.cores = 4L)


names(trees) <- tn

# Building a data set to exctract the information from uniprot
genes <- data.frame(
  tree = rep(tn, sapply(trees, Ntip)),
  gene = unlist(sapply(trees, "[[", "tip.label"))
)

write.csv2(genes, file = "sifter/hundredfamilies_genes.csv", row.names = FALSE)
