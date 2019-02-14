models <- c(
  "01-gold-standard" = "01-gold-standard",
  "02-missing"       = "02-missing",
  "03a-pub-bias"      = "04-full-model",
  "03b-pub-bias-misspecified" = "03-pub-bias"
  )


for (m in models) {
  
  # Listing files
  fns <- list.files(
    path       =  sprintf("simulations/%s", m),
    pattern    =".pdf",
    full.names = TRUE
    )
  
  # New names
  fns_new <- gsub(
    sprintf("simulations/%s/", m),
    sprintf("fig/%s-", m),
    fns,
    fixed = TRUE)
  
  file.copy(from = fns, to = fns_new, overwrite = TRUE)
  
}

# compressing files
zip(
  "fig/fig.zip",
  files = list.files("fig", pattern = ".pdf", full.names = TRUE),
  flags = "-9mj"
  )
