read_tree <- function(file, text) {
  # Removing the end
  x <- if (missing(text) & !missing(file)) readLines(file)
  else if (!missing(text) & missing(file)) text
  else if (!missing(text) & !missing(file)) stop("Choose either of `text` or `file`", call. = FALSE)
  else stop("Either should be specified (`text` or `file`)", call. = FALSE)
  
  # Replacing annotations
  x <- gsub(
    pattern = "(?<=[),(])([a-zA-Z0-9_\\.]*)[:]([0-9\\.]+)\\[&&NHX[:]([a-zA-Z0-9=]+)\\]",
    replacement = "\\1$$$\\3:\\2",
    x = x, perl = TRUE
    )
  
  tree <- ape::read.tree(text=x)
  
  # Annotations
  tip.annotation <- gsub(".*[$]{3}", "", tree$tip.label)
  node.annotation <- gsub(".*[$]{3}", "", tree$node.label)
  tree$tip.label <- gsub("[$]{3}.*", "", tree$tip.label)
  tree$node.label <- gsub("[$]{3}.*", "", tree$node.label)
  
  
  list(tree = tree, tipa = tip.annotation, nodea = node.annotation)
}

# "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F"
# z <- parallel::mclapply(
#   fn_nhx <- list.files("sifter/hundred/", pattern = "nhx", full.names = TRUE),
#   read_tree
# )
# 
# # Counting duplication nodes
# table(unlist(sapply(z, "[[", "nodea")))
# 
# ans <- read_tree("sifter/sulfotransferase/reconciled-pf00685.nhx")
# table(ans$nodea)
# plot(ans$tree, cex=.15)


tree_deminase <- read_tree("sifter/deaminase/reconciled-pf00962-paup.nex")
table(tree_deminase$nodea)

data_deminase <- xml2::read_xml("sifter/deaminase/proteinfamily_pf00962n.pli")
data_deminase <- xml2::as_list(data_deminase)

data_deminase_protein <- data_deminase$Family[names(data_deminase$Family) == "Protein"]
data_deminase_protein <- lapply(data_deminase_protein, data.table::as.data.table)
data_deminase_protein <- data.table::rbindlist(data_deminase_protein, fill = TRUE)
