rm(list = ls())

options(stringsAsFactors = FALSE)
library(aphylo)
library(dplyr)

trees <- list.files("../paint_curated/", full.names = TRUE)

# Checking out an example of data
files <- list.files(trees[2], full.names=TRUE)
names(files) <- gsub(".+[.]", "", files)

#' Function to read GAF files
#' 
#' GO Annotation File Format 2.0
#' 
#' @param x Character scalar. File path
#' @param ... Currently ignored.
#' @return A data frame with 17 columns
#' @references http://www.geneontology.org/page/go-annotation-file-format-20
#' @aliases gaf2.0
read.gaf <- function(x, ...) {
  
  # Reading the file
  dat <- readLines(x)
  
  # EOL should be a space, otherwise R will break it
  dat <- gsub("\t$", "\t ", dat)
  dat <- lapply(dat, strsplit, split="\t", fixed=TRUE)
  
  # Checking version
  if (!grepl("gaf[-]version[:]\\s*2.0", dat[[1]]))
    stop("Currently only GAF version 2.0 is supported.")
  
  dat <- lapply(dat, "[[", 1)[-1]
  
  # Checking lengths
  test <- which(sapply(dat, length) != 17)
  if (length(test))
    stop("Rows ", paste0(test, collapse=", "), " have different length than 17.")
  
  dat <- do.call(rbind, dat)
  colnames(dat) <- c(
    "DB",
    "DB_Object_ID",
    "DB_Object_Symbol",
    "Qualifier",
    "GO_ID",
    "DB:Reference",
    "Evidence_Code",
    "Optional_0",
    "Aspect",
    "DB_Object_Name",
    "DB_Object_Synonym",
    "DB_Object_Type",
    "Taxon",
    "Date",
    "Assigned_By",
    "Annotation_Extension",
    "Gene_Product_Form_ID"
  )
  
  as.data.frame(dat, stringsAsFactors=FALSE)
    
}

# Testing it
pthr_curator  <- read.gaf(files["gaf"])
pthr_attr     <- read.table(files["attr"], sep="\t", header = TRUE)
pthr_exp      <- read.gaf(files["exp"])
pthr_tree     <- read.panther(files["tree"])

# Mapping is:
# Exp(DB + DB_Object_ID) -> Attr(Protein.Id[3]):Attr(Protein.Id[2]) -> tree(tip.label)
mergeit <- function(exp, attr, tree) {
  # Experimental data
  exp_ids  <- with(exp, paste(DB, DB_Object_ID, sep="="))
  
  # Attribute data
  attr_ids <- sapply(sapply(attr[["Protein.Id"]], strsplit, split="\\|"), tail, 1)
  
  # Tree data
  tree_ids <- gsub(".+[:]","",tree$tree$tip.label)
  
  lapply(list(exp=exp_ids, attr=attr_ids, tree=tree_ids), unname)
}
  
ans <- mergeit(pthr_exp, pthr_attr, pthr_tree)