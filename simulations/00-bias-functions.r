bias_calc <- function(fn, dat) {
  
  estimates <- readRDS(fn)
  
  # Extracting the main components, this object will be used throughout the
  # function. At the very end it will be reconciled with the rest of the by-prod
  meta <- estimates %>%
    tibble(tree = .) %>%
    mutate(
      index    = 1:n(),
      is_error = !sapply(tree, inherits, "aphylo_estimates")
    ) %>%
    filter(!is_error) %>%
    mutate(
      estimates = parallel::mclapply(tree, coef, mc.cores=1L),
      variances = parallel::mclapply(tree, function(i) diag(vcov(i)), mc.cores=1L),
      quantiles = parallel::mclapply(tree, function(t_) {
        apply(do.call(rbind, t_$hist), 2L, quantile, probs = c(0.025, .975))
      }, mc.cores=1L),
      NLeafs    = unlist(parallel::mclapply(tree, Ntip)),
      TreeSize  = NLeafs + unlist(parallel::mclapply(tree, Nnode)),
      Missing   = unlist(parallel::mclapply(tree, function(i) sum(i$dat$tip.annotation == 9L), mc.cores=1L)),
      PropOf0   = unlist(parallel::mclapply(tree, function(i) sum(i$dat$tip.annotation == 0L), mc.cores=1L)),
      PropOf0   = PropOf0/(NLeafs - Missing),
      Missing   = Missing/NLeafs
    ) 
  
  # Merging with population data
  coefs_pop <- tibble(
    index = 1L:length(dat),
    par   = lapply(dat, "[[", "par")
  ) %>%
    mutate(is_error = sapply(par, length) == 0) %>%
    filter(!is_error) %>%
    select(-is_error) %$%
    cbind(index, do.call(rbind, par)) %>%
    as_tibble %>%
    mutate(index = as.integer(index))
  
  colnames(coefs_pop)[c(-1, -ncol(coefs_pop))] <- paste0(
    colnames(coefs_pop)[c(-1, -ncol(coefs_pop))], "_pop"
  )
  
  # Extracting coefficients
  coefs_est <- meta %$%
    cbind(index = index, do.call(rbind, estimates)) %>%
    as_tibble %>%
    set_colnames(c("index", paste0(colnames(.)[-1], "_estimated")))
  
  # Extracting variances
  var_est <-  meta %$%
    cbind(index = index, do.call(rbind, variances)) %>%
    as_tibble %>%
    set_colnames(c("index", paste0(colnames(.)[-1], "_var")))
  
  # Extracting lower and upper bounds
  bounds <- meta %>%
    select(index, quantiles) %>%
    mutate(
      lb = lapply(quantiles, "[", i=1, j=),
      ub = lapply(quantiles, "[", i=2, j=)
    ) %>%
    filter(sapply(quantiles, inherits, what="matrix")) %$%
    cbind(
      index = index,
      {do.call(rbind, lb) %>% set_colnames(paste0(colnames(.), "_lb"))},
      {do.call(rbind, ub) %>% set_colnames(paste0(colnames(.), "_ub"))}
    ) %>%
    as_tibble
  
  # Final dataset
  dat_ <- meta %>%
    select(-estimates, -tree, -is_error, -variances, -quantiles) %>%
    left_join(coefs_pop, by = "index") %>%
    left_join(coefs_est, by = "index") %>%
    left_join(var_est, by = "index") %>%
    left_join(bounds, by = "index") %>%
    mutate(index = as.integer(index))
  
  coefs_names <- colnames(dat_)[grepl("[_]estimated$",colnames(dat_))] %>%
    gsub("[_]estimated$", "", .) %>%
    unique()
  
  # Computing coverage
  for (coef_ in coefs_names) {
    
    # Bounds
    dat_[[paste0(coef_, "_covered95")]] <- 
      (dat_[[paste0(coef_, "_lb")]] <= dat_[[paste0(coef_, "_pop")]]) &
      (dat_[[paste0(coef_, "_ub")]] >= dat_[[paste0(coef_, "_pop")]])
    
    dat_[[paste0(coef_, "_bias")]] <-
      dat_[[paste0(coef_, "_estimated")]] - dat_[[paste0(coef_, "_pop")]]
  }

  # Tree size
  dat_$size_tag <- interval_tags(dat_$NLeafs, quantile(dat_$NLeafs, na.rm = TRUE),
                               digits = 0L)

  # NLeafs/TreeSize
  dat_$PropLeafs <- with(dat_, NLeafs/TreeSize)
  dat_$PropLeafs_tag <- interval_tags(dat_$PropLeafs, quantile(dat_$PropLeafs, na.rm=TRUE))

  
  dat_
}

# Creates nice interval tags in the form of [a,b)...[a,z] (last one closed).
# All numbers must be within the interval
interval_tags <- function(x, marks, digits = 1L) {
  
  # Generating labels
  n <- length(marks)
  l <- c(sprintf(paste0("[%.", digits, "f, %.", digits, "f)"), marks[-n][-(n-1)], marks[-n][-1]),
         sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), marks[n-1], marks[n])
  )
  
  # Finding intervals
  x <- findInterval(x, marks, rightmost.closed = TRUE)
  factor(x, levels = 1:length(l), labels = l)
  
}
