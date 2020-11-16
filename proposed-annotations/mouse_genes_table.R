library(data.table)

mouse_genes <- fread("proposed-annotations/mouse_genes.csv")

# Loading proposed predictions to be merged
proposed <- fread("proposed-annotations/proposed-annotations.csv")
proposed[, UniProtKB := gsub(".+=", "", id)]
pretty_table <- merge(
  x = mouse_genes,
  y = proposed, by = "UniProtKB", all.x = TRUE, all.y = FALSE
)

# Filtering data
setnames(
  pretty_table,
  c("Gene Name Gene Symbol", "PANTHER Family/Subfamily"),
  c("Gene Symbol", "PANTHER Family")
  )
pretty_table[, `Gene Symbol` := gsub("^[^\n]+\n|\n[^\n]+$", "", `Gene Symbol`)]
pretty_table[, `PANTHER Family` := gsub(".+\\(|[:]SF[0-9]+\\).*", "", `PANTHER Family`, perl = TRUE)]
pretty_table[, `95\\% C.I.` := sprintf("[%.2f %.2f]", p0_025, p0_975)]
pretty_table[, c(4,6:8) := NULL]
pretty_table[, `Proposed Qualifier` := ifelse(qualifier_proposed == "YES", "", "NOT")]

goterms <- c(
  "GO:0004017" = "adenylate kinase activity",
  "GO:0043065" = "positive regulation of apoptotic process",
  "GO:0004303" = "estradiol 17-beta-dehydrogenase activity",
  "GO:0005615" = "extracellular space",
  "GO:0003677" = "DNA binding",
  "GO:0010039" = "response to iron ion",
  "GO:0004396" = "hexokinase activity",
  "GO:0016298" = "lipase activity",
  "GO:0008543" = "fibroblast growth factor receptor signaling pathway",
  "GO:0004867" = "serine-type endopeptidase inhibitor activity",
  "GO:0008131" = "primary amine oxidase activity",
  "GO:0006888" = "endoplasmic reticulum to Golgi vesicle-mediated transport",
  "GO:0004713" = "protein tyrosine kinase activity",
  "GO:0005739" = "mitochondrion"
  
)

pretty_table[, Annotation := sprintf("%s (%s)",goterms[go_id], go_id)]

pretty_table[, c("go_id", "qualifier_proposed", "PANTHER Family") := NULL]

pretty_table[, Evidence:= ""]

pretty_table
#'     Gene Symbol Qualifier Go_Evidence_Code  95\\% C.I.                                                             Annotation Evidence
#'  1:         Bok                            [0.93 0.99]                  positive regulation of apoptotic process (GO:0043065)         
#'  2:    Serpini1                            [0.94 0.99]              serine-type endopeptidase inhibitor activity (GO:0004867)         
#'  3:        Aoc3                            [0.94 0.99]                            primary amine oxidase activity (GO:0008131)         
#'  4:     Hsd17b1                        IDA [0.93 0.98]                  estradiol 17-beta-dehydrogenase activity (GO:0004303)         
#'  5:       Rab1A                            [0.90 0.98] endoplasmic reticulum to Golgi vesicle-mediated transport (GO:0006888)         
#'  6:        Nfib       NOT                  [0.02 0.10]                                               DNA binding (GO:0003677)         
#'  7:    Pnliprp1       NOT                  [0.02 0.08]                                       extracellular space (GO:0005615)         
#'  8:    Pnliprp1       NOT                  [0.03 0.09]                                           lipase activity (GO:0016298)         
#'  9:    Serpinf2                            [0.93 0.98]              serine-type endopeptidase inhibitor activity (GO:0004867)         
#' 10:       Ephb1                            [0.90 0.97]                                                        NA (GO:0004713)         
#' 11:       Pex14       NOT                  [0.03 0.10]                                                        NA (GO:0007031)         
#' 12:         Ak2                            [0.90 0.97]                                 adenylate kinase activity (GO:0004017)         
#' 13:         Ak3       NOT                  [0.02 0.09]                                 adenylate kinase activity (GO:0004017)         


print(
  xtable::xtable(pretty_table),
  booktabs = TRUE,
  include.rownames = FALSE
)

# Analyzing NOT annotation -----------------------------------------------------
estimates <- readRDS("parameter-estimates/mcmc_partially_annotated_no_prior.rds")
trees     <- readRDS("data/candidate_trees.rds")

# Looking at panther family Human gene O00712 which has a zero
p_with_zero <- predict(estimates, newdata = trees$PTHR11492)
tree_alt <- trees$PTHR11492
tree_alt$tip.annotation[tree_alt$tree$tip.label == "UniProtKB=O00712"] # Is a zero!
tree_alt$tip.annotation[tree_alt$tree$tip.label == "UniProtKB=O00712"] <- 9L

p_without_zero <- predict(estimates, newdata = tree_alt)

# Looking at the mouse gene
p_with_zero["UniProtKB=P97863",]
p_without_zero["UniProtKB=P97863",]
