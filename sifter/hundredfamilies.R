

library(ape)
library(aphylo)

# Reading trees
fn <- list.files("sifter/hundred/", pattern = "nhx", full.names = TRUE)
tn <- gsub("(.+/)([a-zA-Z0-9]+)\\.nhx", "\\2", fn, perl = TRUE)

read_nhx <- function(text) {
  text <- readLines(text)
  
  # This pattern catches most of the information
  pattern <- "([(,)])([a-zA-Z0-9_]*)([:][0-9]+[.]?[0-9]*)(\\[&&NHX[a-zA-Z0-9_:=]+\\])?"
  
  # Capturing the patterns and splitting the data
  x       <- gregexpr(text, pattern = pattern, perl = TRUE)
  x       <- regmatches(text, x)
  
  # Creating a matrix with the data
  x <- regmatches(x[[1]], regexec(x[[1]], pattern = pattern, perl = TRUE))
  x <- do.call(rbind, x)
  
  # Do all have ids?
  noid <- which(x[,3] == "")
  if (length(noid)) {
    x[noid,3] <- sprintf("unnamed%04i", 1:length(noid))
  }
  
  # Replacing the ones with empty id
  for (i in noid) {
    text <- sub(
      pattern     = x[i,1],
      replacement = paste0(x[i,2:4], collapse = ""),
      x           = text,
      fixed       = TRUE
      )
  }
  
  # Is there any root?
  text <- sub(
    pattern = "[)][;]$", replacement = ")root;", x = text, perl = TRUE
    )
  
  dat <- x[,-c(1L, 2L)]
  
  # Capturing NHX fields

  list(
    tree = ape::read.tree(text = text),
    edge = dat[,-3L],
    nhx  = strsplit(dat[,3L], split=":")
  )
  
}

trees <- parallel::mclapply(fn, read_nhx, mc.cores = 4L)
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
