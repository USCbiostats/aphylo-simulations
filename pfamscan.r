#' Process pfamscan result 
#' @param x A string as that returned by EBI pfamscan's API.
#' @return A data.frame with the following columns:
#' - `seq_id`
#' - `alignment_start`
#' - `alignment_end`
#' - `envelope_start`
#' - `envelope_end`
#' - `hmm_acc`
#' - `hmm_name`
#' - `type`
#' - `hmm_start`
#' - `hmm_end`
#' - `hmm_length`
#' - `bit_score`
#' - `E.value`
#' - `significance`
#' - `clan`
pfamscan_process <- function(x) {
  
  x <- strsplit(x, "\\n")[[1]]
  
  # Processing the header
  headers <- x[which(grepl("^\\s*$", x)) - 1]
  headers <- strsplit(headers, ">\\s*<")[[1]]
  headers <- gsub("^#\\s*<|>\\s*$", "", headers)
  headers <- gsub("\\s+", "_", headers)
  
  x <- x[!grepl("^(#.+|\\s*)$", x)]
  
  tmp <- tempfile(fileext = ".txt")
  cat(x, file=tmp, sep = "\n", append=TRUE)
  read.table(
    tmp, header = FALSE, col.names = headers,
    stringsAsFactors = FALSE
  )
}

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
  maxchecktime = 120,
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
    
    ans <- lapply(
      X             = sizes,
      FUN           = pfamscan,
      email         = email,
      ...,
      path          = path,
      maxchecktime  = maxchecktime,
      wait          = wait,
      errorfun      = errorfun,
      max_post_size = max_post_size
    )
    
    return(do.call(rbind, ans))
    
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
  pfamscan_process(httr::content(query3))
  
}