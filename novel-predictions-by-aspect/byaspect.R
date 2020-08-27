library(aphylo)

trees  <- readRDS("data/candidate_trees_inferred.rds")
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

# Did we got the data from QuickGO already?
if (file.exists("novel-predictions-by-aspect/go_terms_info.rds")) {
  
  go_terms_info <- readRDS("novel-predictions-by-aspect/go_terms_info.rds")
  
} else {
  
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
  
  saveRDS(go_terms_info, "novel-predictions-by-aspect/go_terms_info.rds")
  
}

classes <- match(go_terms, go_terms_info$id)
classes <- go_terms_info$aspect[classes]

# Subset of annotations for biological process
write.csv(data.frame(
  go=sapply(trees[classes == "biological_process"], function(i) colnames(i$tip.annotation)),
  tree=tnames[classes == "biological_process"],
  row.names = NULL
), file = "novel-predictions-by-aspect/bio-process.csv", row.names = FALSE, quote = FALSE)

# Fitting the models -----------------------------------------------------------

mle_estimates <- vector("list", 3L)
names(mle_estimates) <- unique(classes)

mcmc_estimates <- mle_estimates
set.seed(123)
for (i in names(mle_estimates)) {
  
  # Finding joint estimates
  tmp_data <- lapply(trees[classes == i], rdrop_annotations, pcent=.85, informative = TRUE)
  tmp_data <- do.call(c, tmp_data)
  mle_estimates[[i]] <- aphylo_mle(tmp_data ~ psi + mu_d + mu_s + Pi)
  
  message("MLE estimates for ", i, " done.")
  
  mcmc_estimates[[i]] <- aphylo_mcmc(
    tmp_data ~ psi + mu_d + mu_s + Pi,
    control = list(
      kernel  = fmcmc::kernel_adapt(lb = .0001, ub = .9999, warmup = 1000, eps = 1e-3),
      nsteps  = 1e5L,
      nchains = 2L,
      thin    = 1L,
      burnin  = 0L,
      conv_checker = NULL,
      progress = TRUE
    )
  )
  
  print(window(mcmc_estimates[[i]], start = 50000))
  
  message("MCMC estimates for ", i, " done.")
  
}

saveRDS(mle_estimates, "novel-predictions-by-aspect/byaspect_mle.rds")
saveRDS(mcmc_estimates, "novel-predictions-by-aspect/byaspect_mcmc.rds")

mle_estimates_ram <- vector("list", 3L)
names(mle_estimates_ram) <- unique(classes)

mcmc_estimates_ram <- mle_estimates_ram
set.seed(123)
for (i in names(mle_estimates_ram)) {
  
  # Finding joint estimates
  tmp_data <- lapply(trees[classes == i], rdrop_annotations, pcent=.85, informative = TRUE)
  tmp_data <- do.call(c, tmp_data)
  mle_estimates_ram[[i]] <- aphylo_mle(tmp_data ~ psi + mu_d + mu_s + Pi)
  
  message("MLE estimates for ", i, " done.")
  
  mcmc_estimates_ram[[i]] <- aphylo_mcmc(
    tmp_data ~ psi + mu_d + mu_s + Pi,
    control = list(
      kernel  = fmcmc::kernel_ram(lb = .0001, ub = .9999, warmup = 1000, eps = 1e-3),
      nsteps  = 1e5L,
      nchains = 2L,
      thin    = 1L,
      burnin  = 0L,
      conv_checker = NULL,
      progress = TRUE
    )
  )
  
  print(window(mcmc_estimates_ram[[i]], start = 50000))
  
  message("MCMC estimates for ", i, " done.")
  
}

saveRDS(mle_estimates_ram, "novel-predictions-by-aspect/byaspect_mle_ram.rds")
saveRDS(mcmc_estimates_ram, "novel-predictions-by-aspect/byaspect_mcmc_ram.rds")

