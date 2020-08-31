library(data.table)

# Reading annotations
goa <- fread("sifter/gene_association.goa_uniprot.80", header = FALSE)

# Setting colnames
cnames <- readLines("sifter/gaf_1.0_colnames.txt")
colnames(goa) <- cnames

# Reading hundred fams
fams <- fread("sifter/hundredfamilies_annotations.csv.gz")

dat <- merge(x = fams, y = goa, by.y = "DB_Object_ID", by.x = "number",
  all.x = TRUE, all.y = FALSE)

saveRDS(dat, "sifter/goa_uniprot80.rds", compress = TRUE)

