library(aphylo)
library(data.table)

# Loading the trees for the experiment
candidate_trees <- readRDS("data/candidate_trees.rds")

# Extracting database with annotations per tree
dat <- lapply(candidate_trees, function(tree) {
  idx <- which(tree$tip.annotation != 9)
  data.frame(
    go        = colnames(tree$tip.annotation),
    UniProtKB = gsub("UniProtKB=", "", tree$tree$tip.label[idx]),
    qualifier = tree$tip.annotation[idx]
  )
})

dat <- Map(
  function(d, n) cbind(d, tree = n),
  d = dat, n = names(candidate_trees)
  )

dat <- rbindlist(dat)

# Obtaining the Pfam name ------------------------------------------------------
# https://www.uniprot.org/help/api_retrieve_entries
get_info_uniprot <- function(uniprotkb, skip_failed = TRUE, db = "pfam") {
  
  # Can be applied to multiple data
  if (length(uniprotkb) > 1L) {
    
    message("Starting multiple queries...")
    
    ans <- vector("list", length(uniprotkb))
    for (i in seq_along(uniprotkb)) {
      
      # Making the query
      ans[[i]] <- get_info_uniprot(
        uniprotkb = uniprotkb[i], skip_failed = skip_failed, db = db
        )
      
      if (!(i %% 10))
        message(sprintf("% 4i/% 4i complete...", i, length(uniprotkb)))
    }
    message("done.")
    names(ans) <- uniprotkb
    return(ans)
  }
  
  # Selecting the source
  query <- if (db == "pfam")
    "http://pfam.xfam.org/protein/%s?output=xml"
  else if (db == "uniprot")
    "https://www.uniprot.org/uniprot/%s.txt"
  
  res <- httr::GET(
    sprintf(query, uniprotkb)
    )
  
  # Checking status
  if (httr::status_code(res) != 200) {
    if (skip_failed) {
      warning("The query failed with status: ", httr::status_code(res))
      return(NULL)
    } else {
      stop("The query failed with status: ", httr::status_code(res))
    }
    
  }
  
  # Returning raw content
  ans <- xml2::as_list(httr::content(res))
  
  if (db == "pfam") {
    return(ans)
  }
  
  return(ans)
  
}

pfam_ids <- unique(dat[,as.character(UniProtKB)])
pfam_ids <- get_info_uniprot(pfam_ids)

save.image("data/sifter_data.rda")

# Setting up the data to be used with SIFTER -----------------------------------

ans <- unique(dat[, .(tree, UniProtKB)])[, .(n=.N), by = tree]
ans[order(n),] # Selecting one tree with small number of proteins annotated

# Extracting the genes (with full names) associated the tree
# PTHR10082
idx <- dat[, which(tree == "PTHR11926")]

sample_trees <- pfam_ids[idx]

cat(names(sample_trees), sep = "\n")

# Creating the fasta sequences for querying pfamscan ---------------------------
fasta <- lapply(sample_trees, function(i) {
  d <- strsplit(
    i$pfam$entry$sequence[[1]],
    split = "(?<=.{50})",
    perl = TRUE)[[1]]
  paste0(
    sprintf(">%s\n", attr(i$pfam$entry, "accession")),
    paste(d, collapse="\n"),
    "*")
})


#' Query the online service PfamScan
#' This uses the EBI API. Details here:
#' https://www.ebi.ac.uk/Tools/common/tools/help/index.html
#' @param fasta_str A string representing a single sequence in fasta format
#' @param email Required by the API
pfamscan <- function(
  fasta_str,
  email,
  ...,
  api_url      = "https://www.ebi.ac.uk/",
  path         = "Tools/services/rest/pfamscan",
  maxchecktime = 60,
  wait         = 1,
  errorfun     = message,
  max_post_size = 2000L
  ) {
  
  # Checking lengths
  if (length(fasta_str) > 1) {
    
    # Learning the size of the sequences
    sizes <- cumsum(nchar(fasta_str)) %/% (max_post_size + 1)
    
    # Making the splits
    sizes <- split(fasta_str, sizes)
    sizes <- sapply(sizes, paste0, collapse = "\n")
    
    return(
      lapply(
        X   = sizes,
        FUN = pfamscan,
        email = email,
        ...,
        path  = path,
        maxchecktime = maxchecktime,
        wait = wait,
        errorfun = errorfun,
        max_post_size = max_post_size
        ))
    
  }

  # Posting the data to the EBI server
  sequences <- gsub("\\n[[:upper:]\\n*]+", "", fasta_str)
  sequences <- strsplit(sequences, "\\n?>")[[1L]][-1]
  message(
    paste(rep("-", 80), collapse = ""),
    sprintf(
      "\nPosting the query for %i sequences:\n %s...",
      length(sequences),
      paste0(sequences, collapse = ", ")
      ),
    appendLF = TRUE
    )
  query1 <- httr::POST(
    url    = api_url,
    path   = c("Tools/services/rest/pfamscan", "run"),
    config = httr::add_headers(
      "Content-Type" = "application/x-www-form-urlencoded",
      Accept         = "text/plain"
      ),
    encode ="form",
    body = list(
      email    = email,
      database = "pfam-a",
      format   = "txt",
      sequence = fasta_str
      ),
    ...
  )
  
  # Checking status
  if (httr::status_code(query1) != 200) {
    errorfun("Query did not worked: ", httr::content(query1))
    return(NULL)
  } 
  
  # Retrieving the id
  query_id <- httr::content(query1)
  
  # Checking if it is done
  time0 <- Sys.time()
  message("Waiting for the query: ", query_id, " to finalize...", appendLF = FALSE)
  while (difftime(Sys.time(), time0, units="secs") < maxchecktime) {
    
    query2 <- httr::GET(
      url = api_url,
      path = c("Tools/services/rest/pfamscan", "status", query_id),
      httr::add_headers(
        Accept = "text/plain"
      ),
      ...
    )
    
    # Checking the status
    if (httr::status_code(query2) != 200) {
      errorfun(
        sprintf(
          "Querying the id %s did not worked: %s",
          query_id,
          httr::content(query2)
        )
      )
      
      return(NULL)
    }
    
    # Checking whether it is done or not
    if (httr::content(query2) == "FINISHED")
      break
    
    message(".", appendLF = FALSE)
    
    Sys.sleep(wait)
    
  }
  
  if (difftime(Sys.time(), time0, units="secs") >= maxchecktime) {
    errorfun("Checking query ", query_id, " went out-of-time.")
    return(NULL)
  }
    
  
  message("Success!")
  
  # Retrieving the results
  # curl -X GET --header 'Accept: text/plain' ''
  message("Retrieving the results for query ", query_id, "...", appendLF = )
  query3 <- httr::GET(
    url = api_url,
    path = c("Tools/services/rest/pfamscan", "result", query_id, "out"),
    httr::add_headers(
      Accept = "text/plain"
    ),
    ...
  )
  
  if (httr::status_code(query3) != 200) {
    errorfun(
      sprintf(
        "Retrieving the query id %s did not worked: %s",
        query_id,
        httr::content(query3)
      )
    )
    return(NULL)
  }
  message("Done!")
  httr::content(query3)

}

# pfamscan_results <- vector("list", length(fasta))
# for (i in 1:length(pfamscan_results)) {
#   pfamscan_results[[i]] <- pfamscan(
#     fasta_str = fasta[[i]],
#     email     = "vegayon@usc.edu"
#     )
# }
# 
# pfamscan_results_cleaned <- unlist(pfamscan_results)
# pfamscan_results_cleaned <- unlist(strsplit(pfamscan_results_cleaned, split="\n"))
# pfamscan_results_cleaned <- pfamscan_results_cleaned[
#   grepl("^[^#]", pfamscan_results_cleaned)
#   ]

pfamscan_results2 <- pfamscan(
  fasta_str     = fasta,
  email         = "vegayon@usc.edu",
  max_post_size = 50000
)

saveRDS(pfamscan_results, "data/sifter_data_pfamscan.rds")
# 
# 
# # Preparing the command
# sprintf(
#   "python sifter_find_families.py -p %s %s",
#   paste(names(sample_trees),collapse=","),
#   "../aphylo_family_list.txt"
#   )
# 
# sprintf(
#   "python sifter_gather_family_data.py -i %s %s" ,
#   "../aphylo/family_list.txt",
#   "../aphylo/fam_data"
# )
# 
# 
# 
# sample_trees$Q07813$pfam$entry

