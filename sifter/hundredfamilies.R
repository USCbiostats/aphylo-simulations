

library(ape)
library(aphylo)
library(stringi)

# Reading trees
fn <- list.files("sifter/hundred/", pattern = "nhx", full.names = TRUE)
tn <- gsub("(.+/)([a-zA-Z0-9]+)\\.nhx", "\\2", fn, perl = TRUE)

read_nhx <- function(text) {
  text <- readLines(text)
  
  # This pattern catches most of the information
  pattern <- "([(,)])([a-zA-Z0-9_]*)([:][0-9]+[.]?[0-9]*)(\\[&&NHX[a-zA-Z0-9_:=]+\\])?"
  
  # Capturing the patterns and splitting the data
  x <- gregexpr(text, pattern = pattern, perl = TRUE)
  x <- regmatches(text, x)
  
  # Creating a matrix with the data
  x <- regmatches(x[[1]], regexec(x[[1]], pattern = pattern, perl = TRUE))
  x <- do.call(rbind, x)
  
  # Do all have ids?
  noid <- which(x[,3] == "")
  if (length(noid)) {
    x[noid,3] <- sprintf("unnamed%04i", 1:length(noid))
  }
  
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
  nhx <- strsplit(dat[,3], split = "[:]|[=]")
  nhx <- tryCatch(lapply(nhx, function(n) {
    if (length(n) > 0) {
      n <- gsub(pattern = "(^\\[|\\]$)", replacement = "", x = n)
      n <- matrix(n[-1], ncol = 2, byrow = TRUE)
      structure(.Data = n[,2], names = n[,1])
    } else
     n
  }), error = function(e) e)
  
  if (inherits(nhx, "error")) 
    stop(
      "There was a problem when processing the &&NHS blocks.",
      " Possibly, not all the attributes have the right tag. Here is the error",
      ":\n", paste0(nhx, collapse=""), call. = FALSE
      )
  
  list(
    tree = ape::read.tree(text = text),
    edge = dat[,-3L],
    nhx  = nhx
  )
  
}

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
