

library(ape)
library(aphylo)

# Reading trees
fn_nhx <- list.files("sifter/hundred/", pattern = "nhx", full.names = TRUE)
tn <- gsub("(.+/)([a-zA-Z0-9]+)\\.nhx", "\\2", fn, perl = TRUE)

trees <- parallel::mclapply(fn, read.tree, mc.cores = 4L)
info  <- parallel::mclapply(fn, function(i) {
  x <- readLines(i)
  stringr::str_extract_all(x, "\\[:\\][0-9\\.]+\\[&&NHX:[A-Z][:=][A-Z]\\]")
}, mc.cores = 4L)

names(trees) <- tn

# Reading functional annotations
fn_pli