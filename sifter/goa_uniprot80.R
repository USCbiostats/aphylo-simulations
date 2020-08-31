library(data.table)

# Reading annotations
goa <- fread("sifter/gene_association.goa_uniprot.80", header = FALSE)

# Setting colnames
cnames <- readLines("sifter/gaf_1.0_colnames.txt")
colnames(goa) <- cnames

# Keeping only relevant
codes <- c("IDA", "IMP", "TAS")
goa <- goa[Evidence_Code %in% codes]

# ID Synonim
goa[, synm := gsub(".+[|]([a-zA-Z0-9_]+)$", "\\1", DB_Object_Synonym, perl = TRUE)]

# Reading hundred fams
fams <- fread("sifter/hundredfamilies_annotations.csv")

dat <- merge(x = fams, y = goa, by.y = "synm", by.x = "number",
  all.x = TRUE, all.y = FALSE)

saveRDS(dat, "sifter/goa_uniprot80.rds", compress = TRUE)

