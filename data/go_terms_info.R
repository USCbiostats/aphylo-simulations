library(aphylo)

trees <- readRDS("data/candidate_trees.rds")

tnames <- names(trees)
tnames <- rep(tnames, sapply(trees, Nann))

trees <- do.call(c, trees)
trees <- unlist(lapply(trees, function(d) {
  lapply(1:Nann(d), function(i) d[,i])
}), recursive = FALSE)
trees <- do.call(c, trees)

# Figuring out the aspect of each function
# (biological process, cellular component, molecular function)
go_terms <- sapply(trees, function(i) colnames(i$tip.annotation))

quick_go_call <- function(geneProductId) {
  
  # API call baseline
  baseline_url  <- "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"
  geneProductId <- paste(unique(geneProductId), collapse = "%2C")
  
  query <- sprintf("-X GET --header 'Accept:application/json' '%s%s'", baseline_url, geneProductId)
  message("Submitting query:\n", query)
  
  system2("curl", query, stdout = TRUE)
  
}

go_terms_info <- quick_go_call(unique(go_terms))
go_terms_info <- jsonlite::fromJSON(go_terms_info)
go_terms_info <- go_terms_info$results[
  ,c("id", "name", "definition", "aspect")
]

go_terms_info$definition <- go_terms_info$definition$text


which(!unique(go_terms) %in% go_terms_info$id)
# 114, which is the term "GO:0004882", which is tagged as obsolete.

go_terms_info <- rbind(
  go_terms_info,
  data.frame(
    id         = "GO:0004882",
    name       = "nuclear receptor activity",
    definition = "Combining with a signal and transmitting the signal to the transcriptional machinery by interacting selectively and non-covalently with a specific double-stranded genomic DNA sequence in order to modulate transcription by RNA polymerase II.",
    aspect     = "molecular_function",
    stringsAsFactors = FALSE
  )
)


data.table::fwrite(go_terms_info, "data/go_terms_info.csv")
