library(aphylo)
source("sifter/read_sifter.R")

tree <- read_nhx("sifter/deaminase/reconciled-pf00962-paup.nex")
ann  <- read_pli("sifter/deaminase/proteinfamily_pf00962n.pli")
x <- xml2::read_html("sifter/deaminase/proteinfamily_pf00962n.pli")
x <- xml2::as_list(x)
2