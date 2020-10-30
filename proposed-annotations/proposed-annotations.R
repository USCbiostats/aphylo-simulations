library(aphylo)

estimates <- readRDS("parameter-estimates/mcmc_partially_annotated_no_prior.rds")

# Predictions
divisions <- parallel::splitIndices(Ntrees(estimates), 4)
cl <- parallel::makeForkCluster(4)
ans <- parallel::parLapply(cl, divisions, function(i) {
  predict(
    object     = estimates,
    which.tree = i,
    loo        = FALSE,
    nsamples   = 100, 
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

# Selecting highly accurate predictions
dat <- subset(dat, p0_025 > .90 | p0_975 < .10)

# That are not in the original set
original <- sapply(estimates$dat, function(d) {
  
  d$tree$tip.label[d$tip.annotation != 9]
  
})
original <- unlist(original)

# First subset
dat <- subset(dat, !(id %in% original))

library(data.table)
inferred <- fread("data-raw/inferred_annotations")
inferred[, id := gsub(".+UniProtKB", "UniProtKB", V2)]
inferred <- inferred[, list(id, V3)]
inferred[, exists := TRUE]

discoveries <- merge(x = dat, y = inferred, by.x = c("id", "fun"), by.y = c("id", "V3"), all.x = TRUE, all.y = FALSE)

# Taking one example
which(sapply(estimates$dat, function(i) "UniProtKB=A0A0B4J2U9" %in% i$tree$tip.label))

graphics.off()
pdf("parameter-estimates/proposed-annotations.pdf", width = 9, height = 9)
plot(
  estimates, which.tree = 43,
  nsamples = 500, loo = FALSE,
  ncores = 4, cex = .1
  )
dev.off()

# Expornting table
discoveries <- as.data.table(discoveries)
discoveries <- discoveries[is.na(exists) == TRUE]

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
    c(selectedFields, sprintf("geneProductId=%s", geneProductId)),
    collapse = "&"
  )
  
  query <- sprintf("-X GET --header 'Accept:text/gpad' '%s%s'", baseline_url, query)
  message("Submitting query:\n", query)
  
  system2("curl", query, stdout = TRUE)
  
}

# Getting the first 20
# quick_go <- vector("list", 20)
discoveries[, idsolo := gsub(".+=", "", id)]
quickgo_data <- quick_go_call(discoveries$idsolo)
quickgo_data <- fread(text = quickgo_data, skip = 10)

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

in_go <- quickgo_data[, list(DB_Object_ID, GO_ID)]
in_go[, quick_go := TRUE]

discoveries <- merge(
  discoveries,
  y = in_go,
  by.x = c("id", "fun"),
  by.y = c("DB_Object_ID", "GO_ID"),
  all.x = TRUE,
  all.y = FALSE
)

discoveries <- subset(discoveries, select = c(-exists, -idsolo, -quick_go))
colnames(discoveries)[1:2] <- c("id", "go_id")
discoveries[, annotation := ifelse(p0_025 > .8, "YES", "NOT")]

data.table::fwrite(
  discoveries,
  file = "parameter-estimates/proposed-annotations.csv"
  )
