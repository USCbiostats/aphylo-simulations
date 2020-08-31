

library(ape)
library(aphylo)

# Reading trees
fn <- list.files("sifter/sifter/100trees/", pattern = "nhx", full.names = TRUE)
tn <- gsub("(.+/)([a-zA-Z0-9]+)\\.nhx", "\\2", fn, perl = TRUE)

trees <- parallel::mclapply(fn, read.tree, mc.cores = 4L)
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
