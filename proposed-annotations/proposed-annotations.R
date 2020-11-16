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

fwrite(dat, "proposed-annotations/proposed-annotations-raw.csv.gz", compress = 'gzip')

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
      sprintf("geneProductId=%s", geneProductId) #,
      # "evidenceCode=ECO:0000006,ECO:0006055&evidenceCodeUsage=descendants"
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

in_go <- copy(quickgo_data)

in_go[, go_evidence_code := gsub("goEvidence=([A-Z]+)", "\\1", Annotation_Properties, perl = TRUE)]
in_go[, qualifier_binary := fifelse(grepl("^NOT", Qualifier), 0L, 1L)]

in_go <- in_go[, list(DB_Object_ID, GO_ID, qualifier_binary, go_evidence_code)]

in_go[, go_evidence_code := paste(go_evidence_code, collapse=","), by = .(DB_Object_ID, GO_ID, qualifier_binary)]
in_go <- unique(in_go[, .(DB_Object_ID,GO_ID,qualifier_binary,go_evidence_code)])
discoveries[, UniProtKB := gsub(".+=", "", id)]

after_merge <- merge(
  unique(discoveries),
  y = unique(in_go),
  by.x = c("UniProtKB", "fun"),
  by.y = c("DB_Object_ID", "GO_ID"),
  all.x = TRUE,
  all.y = FALSE
)

after_merge[,c("idsolo", "id", "score") := NULL]
colnames(after_merge)[1:2] <- c("UniProtKB", "go_id")
after_merge[, qualifier_proposed := ifelse(p0_025 > .8, "YES", "NOT")]

after_merge[, table(qualifier_proposed, qualifier_binary )]
#                    Qualifier_Binary
# Qualifier_proposed   0   1
#                NOT  10  15
#                YES   0 154

# How many duplicated (contradicting)
table(after_merge[, table(paste0(UniProtKB, go_id))])
#  1   2 
#215   5

# Which are the ones contradicting
# All the cases contradicting, we propose a negative annotation
after_merge[, opposing := .N > 1, by = .(UniProtKB, go_id)][opposing == TRUE]
#'       UniProtKB      go_id     p0_025     p0_500     p0_975 Qualifier_Binary Go_Evidence_Code Qualifier_proposed opposing
#'  1: A0A0G2K1M7 GO:0005229 0.02406897 0.05078447 0.09844351                0          ISO,ISO                NOT     TRUE
#'  2: A0A0G2K1M7 GO:0005229 0.02406897 0.05078447 0.09844351                1          ISO,ISO                NOT     TRUE
#'  3:     F1LY14 GO:0005229 0.02393068 0.05051724 0.09813917                1          IBA,ISO                NOT     TRUE
#'  4:     F1LY14 GO:0005229 0.02393068 0.05051724 0.09813917                0          ISO,ISO                NOT     TRUE
#'  5:     F1LZ77 GO:0005229 0.02408001 0.05085815 0.09882460                0          ISO,ISO                NOT     TRUE
#'  6:     F1LZ77 GO:0005229 0.02408001 0.05085815 0.09882460                1              ISO                NOT     TRUE
#'  7:     F1M0A0 GO:0005229 0.02392630 0.05062639 0.09856256                1              IBA                NOT     TRUE
#'  8:     F1M0A0 GO:0005229 0.02392630 0.05062639 0.09856256                0          ISO,ISO                NOT     TRUE
#'  9:     Q76HM9 GO:0004415 0.02308549 0.05012509 0.08740170                0          ISO,ISO                NOT     TRUE
#' 10:     Q76HM9 GO:0004415 0.02308549 0.05012509 0.08740170                1      ISO,ISO,ISS                NOT     TRUE

# Experimental codes
experimental_codes <- c(
  "EXP",
  "IDA",
  "IPI",
  "IMP",
  "IGI",
  "IEP",
  "HTP",
  "HDA",
  "HMP",
  "HGI",
  "HEP"
  )
experimental_codes <- paste0(experimental_codes, collapse = "|")

# How much agreement
after_merge[opposing != TRUE, table(qualifier_binary, qualifier_proposed, useNA = "always")]

after_merge[
  opposing != TRUE & grepl(experimental_codes, go_evidence_code)
  ]
#    UniProtKB      go_id     p0_025    p0_500     p0_975 qualifier_binary        go_evidence_code qualifier_proposed opposing
# 1:    P51656 GO:0004303 0.92328093 0.9642971 0.98287054                1 IDA,IEA,IEA,ISO,ISO,IEA                YES    FALSE
# 2:    Q9DC29 GO:0005739 0.02185941 0.0522652 0.09772889                1     IEA,HDA,ISO,IEA,IEA                NOT    FALSE

data.table::fwrite(
  after_merge,
  file = "proposed-annotations/proposed-annotations.csv"
  )

# Plot -------------------------------------------------------------------------
# after_merge <- data.table::fread("proposed-annotations/proposed-annotations.csv")
dat <- copy(after_merge)
dat[, novel := is.na(qualifier_binary)]
dat[, exp_evidence := grepl(experimental_codes, go_evidence_code)]
dat <- unique(dat[, .(UniProtKB, go_id, opposing, exp_evidence, novel)])

dat[, note := fifelse(
  exp_evidence, "EXP evidence",
  fifelse(
    opposing, "Dual annotation (Present/Absent)",
    fifelse(
      novel, "Not in GO",
      "non-EXP evidence"
    )
  )
)]

library(ggplot2)    
toplot <- dat[, .(`N Annotations` = .N), by = note]

graphics.off()
pdf("proposed-annotations/proposed-annotations.pdf", width=5, height = 3)
ggplot(toplot, aes(x = note, y = `N Annotations`, fill=note)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = `N Annotations`, y = `N Annotations`), vjust = -.2
    ) + 
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "Number of Annotations", fill = "Legend") +
  theme_minimal(base_family = "serif") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
dev.off()
