library(aphylo)
library(data.table)

estimates <- readRDS("parameter-estimates/mcmc_partially_annotated_no_prior.rds")

# Predictions
divisions <- parallel::splitIndices(Ntrees(estimates), 4)
set.seed(123)
cl <- parallel::makeForkCluster(10)
ans <- parallel::parLapply(cl, divisions, function(i) {
  predict(
    object     = estimates,
    which.tree = i,
    loo        = FALSE,
    nsamples   = 400, 
    ncores     = 1L
)
})
parallel::stopCluster(cl)

# Combining cases
ans <- do.call(c, ans)
  
# Removing NAs (internal nodes)
ans <- lapply(ans, function(a) a[complete.cases(a),])

# Collapsing dataset
dat <- lapply(ans, function(a) {
  data.frame(
    id   = rownames(a),
    p0_025 = a[,1],
    p0_500 = a[,2],
    p0_975 = a[,3],
    fun    = gsub("_.+", "", colnames(a)[1]),
    stringsAsFactors = FALSE
  )
})

dat <- do.call(rbind, dat)
dat <- data.table(dat)

# Selecting highly accurate predictions
dat <- dat[p0_025 > .9 | p0_975 < .1]

# That are not in the original set
original <- lapply(estimates$dat, function(d) {
  
  data.table(
    id = d$tree$tip.label[d$tip.annotation != 9],
    go = colnames(d$tip.annotation)
  )
  
})
original <- rbindlist(original)
original[, in_original := TRUE]

# First subset
discoveries <- merge(
  dat,
  original,
  by.x = c("id", "fun"),
  by.y = c("id", "go"),
  all.x = TRUE, all.y = FALSE
  )

discoveries <- discoveries[is.na(in_original)][, in_original := NULL]

discoveries[, score := ifelse(p0_025 > .5, p0_025, 1 - p0_975)]
discoveries <- discoveries[order(score, decreasing = TRUE)]


# Curl call
quick_go_call <- function(
  geneProductId,
  selectedFields = NULL
  ) {
  
  # API call baseline
  baseline_url  <- "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?"
  if (!is.null(selectedFields))
    selectedFiels <- sprintf("selectedFields=%s", selectedFields)
  geneProductId <- paste(unique(geneProductId), collapse = "%2C")
  query <- paste(
    c(
      selectedFields,
      sprintf("geneProductId=%s", geneProductId),
      "evidenceCode=ECO:0000006,ECO:0006055&evidenceCodeUsage=descendants"
      ),
    collapse = "&"
  )
  
  query <- sprintf("-X GET --header 'Accept:text/gpad' '%s%s'", baseline_url, query)
  message("Submitting query:\n", query)
  
  system2("curl", query, stdout = TRUE)
  
}

# Getting the first 20
# quick_go <- vector("list", 20)
discoveries[, idsolo := gsub(".+=", "", id)]
quickgo_data <- quick_go_call(unique(discoveries$idsolo))
quickgo_data <- fread(text = quickgo_data[grep("^[^!]", quickgo_data)])

colnames(quickgo_data) <- c(
  "DB",
  "DB_Object_ID",
  "Qualifier",
  "GO_ID",
  "DBReference",
  "Evidence_Code",
  "With_or_From",
  "Interacting_taxon_ID",
  "Date",
  "Assigned_by",
  "Annotation_Extension",
  "Annotation_Properties"
)

quickgo_data[, Go_Evidence_Code := gsub("goEvidence=([A-Z]+)", "\\1", Annotation_Properties, perl = TRUE)]

in_go <- quickgo_data[, list(DB_Object_ID, GO_ID, Qualifier, Go_Evidence_Code)]
in_go[, quick_go := TRUE]

discoveries[, UniProtKB := gsub(".+=", "", id)]

after_merge <- merge(
  discoveries,
  y = in_go,
  by.x = c("UniProtKB", "fun"),
  by.y = c("DB_Object_ID", "GO_ID"),
  all.x = TRUE,
  all.y = FALSE
)

after_merge[, table(fifelse(p0_025 < .5, "NOT", "YES"), Qualifier )]

after_merge <- unique(after_merge)

after_merge[,c("idsolo", "quick_go", "id", "score") := NULL]
colnames(after_merge)[1:2] <- c("id", "go_id")
after_merge[, annotation := ifelse(p0_025 > .8, "YES", "NOT")]

data.table::fwrite(
  after_merge,
  file = "proposed-annotations/proposed-annotations.csv"
  )
