

library(ape)
library(aphylo)

# Reading trees
fn <- list.files("sifter/hundred/", pattern = "nhx", full.names = TRUE)
tn <- gsub("(.+/)([a-zA-Z0-9]+)\\.nhx", "\\2", fn, perl = TRUE)

trees <- parallel::mclapply(fn, read.tree, mc.cores = 4L)
info  <- parallel::mclapply(fn, function(i) {
  x <- readLines(i)
  stringr::str_extract_all(x, "\\[:\\][0-9\\.]+\\[&&NHX:[A-Z][:=][A-Z]\\]")
}, mc.cores = 4L)

pf04988 <- readLines("sifter/hundred/PF04988.nhx")
p1 <- "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1[&&NHX:AA=1:B=2])F;"

p <- p1

pattern <- "([a-zA-Z0-9_]*)([:][0-9]+[.]?[0-9]*)(\\[&&NHX([:][a-zA-Z0-9_]+[=][a-zA-Z0-9_]+)*\\])?"

# https://regex101.com/r/oDp844/1
# ([a-zA-Z0-9_]*)([:][0-9]+[.]?[0-9]*)(\[&&NHX[a-zA-Z0-9_:=]+\])?
x       <- gregexpr(p, pattern = pattern, perl = TRUE)
x       <- regmatches(p, x)

regmatches(x[[1]], regexec(x[[1]], pattern = pattern, perl = TRUE))

read_nhx <- function(x) {
  
}


names(trees) <- tn

# Building a data set to exctract the information from uniprot
genes <- data.frame(
  tree = rep(tn, sapply(trees, Ntip)),
  gene = unlist(sapply(trees, "[[", "tip.label"))
)

write.csv2(genes, file = "sifter/hundredfamilies_genes.csv", row.names = FALSE)
