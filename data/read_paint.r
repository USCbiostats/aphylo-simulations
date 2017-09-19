rm(list = ls())

options(stringsAsFactors = FALSE)
library(aphylo)
library(dplyr)
library(parallel)


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
  
  # Removing comment lines
  dat <- lapply(dat, "[[", 1)
  dat <- dat[!grepl("^[!]", dat)]
  
  if (!length(dat))
    stop("The file -", x, "- only had comments.")
  
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
  
  # I Apply this (unique) since it seems that there are duplicates
  unique(as.data.frame(dat, stringsAsFactors=FALSE))
    
}

# Mapping is:
# Exp(DB + DB_Object_ID) -> Attr(Protein.Id[3]):Attr(Protein.Id[2]) -> tree(tip.label)
mergeit <- function(path) {
  
  # Checking out an example of data
  files <- list.files(path, full.names=TRUE)
  names(files) <- gsub(".+[.]", "", files)
  
  substr("'Hola'", 1, 1)
  
  # Reading the data files
  # pthr_curator  <- read.gaf(files["gaf"])
  attr       <- readLines(files["attr"]) # read.table(files["attr"], sep="\t", header = TRUE, quote="")
  attrcnames <- strsplit(attr[1], split="\t")[[1]]
  attrcnames <- sapply(attrcnames, function(x) substr(x, 2L, nchar(x) - 1L))
  attrcnames <- stringr::str_replace_all(attrcnames, " ", ".")
    
  attr       <- do.call(rbind, lapply(lapply(attr[-1L], strsplit, split = "\t"), "[[", 1L))
  colnames(attr) <- attrcnames
  attr       <- apply(attr, 2, function(x) substr(x, rep(2L, length(x)), nchar(x) - 1L))
  attr       <- as.data.frame(attr, stringsAsFactors=FALSE)

  exp      <- read.gaf(files["exp"])
  tree     <- read.panther(files["tree"])
  
  # Experimental data
  # We are only interested in those that have a qualifier
  exp <- subset(exp, Qualifier != "")
  
  # Creating an id to match with attribute data
  exp <- with(
    exp, 
    data.frame(
      exp_to_attr_id = paste(DB, DB_Object_ID, sep="="),
      qualifier      = Qualifier,
      go_id          = GO_ID 
      )
  )
  
  # Attribute data
  attr_ids <- sapply(attr[["Protein.Id"]], strsplit, split="\\|")
  attr_ids <- unname(do.call(rbind, attr_ids))
  
  attr <- data.frame(
    exp_to_attr_id  = attr_ids[,3],
    attr_to_tree_id = attr_ids[,2]
  )
  
  # Tree data
  tree <- data.frame(
    attr_to_tree_id = gsub(".+[:]","",tree$tree$tip.label)
  )
  
  # Merging data ---------------
  ans <- merge(exp, attr, by = "exp_to_attr_id", all.x=TRUE, all.y=FALSE)
  ans <- merge(tree, ans, by = "attr_to_tree_id", all.x=TRUE, all.y=FALSE)
  
  with(ans, data.frame(id = attr_to_tree_id, annotation = qualifier, go_id= go_id))
}

tryMergeit <- function(...) tryCatch(mergeit(...), error = function(e) e)

# Processing the data ----------------------------------------------------------

trees <- list.files("../paint_curated", pattern = "PTHR[0-9]+",full.names = TRUE)

cl <- makeForkCluster(10)

intervals <- seq(1, length(trees), by = 400)
ans <- NULL
for (i in 2:length(intervals)) {
  ans <- c(ans, parLapply(cl, trees[intervals[i-1]:(intervals[i]-1)], tryMergeit))
  message(sprintf("Interval %i to %i complete", intervals[i-1], intervals[i]-1))
}

ans <- c(ans, parLapply(cl, trees[intervals[i]:length(trees)], tryMergeit))
names(ans) <- trees

# How many errors
ans_class <- unlist(lapply(lapply(ans, class), paste, collapse=" "))
table(ans_class)

# Taking a look at errors
which_error <- which(unlist(sapply(ans, inherits, what="error")))
head(ans[which_error])

stopCluster(cl)

# Analyzing the data -----------------------------------------------------------

# Keeping GO+PANTHER+GENE that have a value
tabs <- lapply(trees, function(x) {
  if (inherits(ans[[x]], "error"))
    return(NULL)
  
  d <- ans[[x]][complete.cases(ans[[x]]),]
  if (!nrow(d))
    return(NULL)
  data.frame(pthr = x, d)
})

# Setting annotations to TRUE/FALSE
tabs <- do.call(rbind, tabs)
tabs[["annotation"]] <- !grepl("^(N|n)", tabs[["annotation"]])

stats <- tabs %>% group_by(pthr, go_id) %>%
  summarize(
    yes = sum(annotation),
    no  = sum(!annotation)
  )

# Which have both TRUE/FALSE
View(stats[stats$yes >0 & stats$no > 0,], "Candidates")
# A tibble: 6 x 4
# Groups:   pthr [5]
#                        pthr      go_id   yes    no
#                        <chr>      <chr> <int> <int>
# 1 ../paint_curated/PTHR11361 GO:0032137     4     1
# 2 ../paint_curated/PTHR11361 GO:0032357     4     1
# 3 ../paint_curated/PTHR15137 GO:0005669     2     2
# 4 ../paint_curated/PTHR18898 GO:0000776     2     2
# 5 ../paint_curated/PTHR24115 GO:0005876     4     4
# 6 ../paint_curated/PTHR24256 GO:0031012     2     2
