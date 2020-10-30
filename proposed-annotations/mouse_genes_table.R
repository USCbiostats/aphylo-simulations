library(data.table)

mouse_genes <- fread("proposed-annotations/mouse_genes.csv")

# Loading proposed predictions to be merged
proposed <- fread("proposed-annotations/proposed-annotations.csv")
proposed[, id := gsub(".+=", "", id)]

pretty_table <- merge(
  x = mouse_genes,
  y = proposed, by.x = "UniProtKB", by.y = "id", all.x = TRUE, all.y = FALSE
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
pretty_table[, `Qualifier` := ifelse(annotation == "YES", "", "NOT")]

goterms <- c(
  "GO:0004017" = "adenylate kinase activity",
  "GO:0043065" = "positive regulation of apoptotic process",
  "GO:0004303" = "estradiol 17-beta-dehydrogenase activity",
  "GO:0005615" = "extracellular space",
  "GO:0003677" = "DNA binding",
  "GO:0010039" = "response to iron ion",
  "GO:0004396" = "hexokinase activity",
  "GO:0016298" = "lipase activity"
)

pretty_table[, Annotation := sprintf("%s (%s)",goterms[go_id], go_id)]

pretty_table[, c("go_id", "UniProtKB", "annotation", "PANTHER Family") := NULL]

pretty_table[, Evidence:= ""]

pretty_table

print(
  xtable::xtable(pretty_table),
  booktabs = TRUE,
  include.rownames = FALSE
)
