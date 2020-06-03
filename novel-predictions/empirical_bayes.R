library(aphylo)

trees <- readRDS("data/candidate_trees.rds")
trees <- do.call(c, trees)
trees <- unlist(lapply(trees, function(d) {
  lapply(1:Nann(d), function(i) d[,i])
}), recursive = FALSE)
trees <- do.call(c, trees)

# Figuring out the aspect of each function
# (biological process, cellular component, molecular function)
go_terms <- sapply(trees, function(i) colnames(i$tip.annotation))

# Did we got the data from QuickGO already?
if (file.exists("novel-predictions/go_terms_info.rds")) {
  
  go_terms_info <- readRDS("novel-predictions/go_terms_info.rds")
  
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
  
  saveRDS(go_terms_info, "novel-predictions/go_terms_info.rds")
  
}

classes <- match(go_terms, go_terms_info$id)
classes <- go_terms_info$aspect[classes]
classes[122] <- "molecular_function" # For some reason is NA

if (anyNA(classes))
  stop("One or more GO-terms has no aspect. Need to download the data again!")


# Estimating empirical bayes using map
ALPHAS <- c(psi0 = 2, psi1 = 2, mu_d0 = 18, mu_d1 = 10, mu_s0 = 2, mu_s1 = 2, Pi = 2)
BETAS  <- c(psi0 = 18, psi1 = 18, mu_d0 = 2, mu_d1 = 10, mu_s0 = 18, mu_s1 = 18, Pi = 18)

bpriors <- bprior(shape1 = ALPHAS, shape2 = BETAS)

# How evenly distributed are the types?
prop_types <- sapply(trees, function(i) mean(with(i, c(tip.type, node.type)[reduced_pseq])))
tapply(X = prop_types, INDEX = classes, mean)


# Fitting models per class of function -----------------------------------------
mle_per_aspect <- vector("list", length(unique(classes)))
names(mle_per_aspect) <- unique(classes)
mcmc_per_aspect <- mle_per_aspect
set.seed(1234)
for (i in names(mcmc_per_aspect)) {
  
  # Obtaininng mles first
  mle_per_aspect[[i]] <- aphylo_mle(
    trees[classes == i] ~ psi + mu_d + mu_s + Pi,
    control = list(factr = 1e3) # Smaller factor means more steps, the default value
                                # is 1e7. This is multiplied by .Machine$double.eps
    )
  
  # Obtaining the MCMC estimates using the MLE's as starting point
  mcmc_per_aspect[[i]] <- aphylo_mcmc(
    trees[classes == i] ~ psi + mu_d + mu_s + Pi,
    params  = coef(mle_per_aspect[[i]]),
    # The RAM kernel is not behaving as expected. AM does better and actually 
    # gets the acceptance rate closer to 23%, at least in the case of 
    # molecular_function.
    control = list(
      kernel = fmcmc::kernel_am(lb = .0001, ub = .9999, warmup = 1000, eps = 1e-3),
      nsteps  = 5e4L,
      nchains = 4L,
      thin    = 1L,
      burnin  = 0L,
      conv_checker = NULL,
      progress = TRUE
      )
    )
  message("Tree ",i, " done.")
}

saveRDS(mle_per_aspect, "novel-predictions/empirical_bayes_mle.rds")
saveRDS(mcmc_per_aspect, "novel-predictions/empirical_bayes_mcmc.rds")

# Biological Process the MCMC is doing way better, so it is a good idea to
# rerun the MLE starting from the MCMC estimates
mle_per_aspect$biological_process <- aphylo_mle(
  trees[classes == "biological_process"] ~ psi + mu_d + mu_s + Pi,
  params = coef(window(mcmc_per_aspect$biological_process, start = 30000)),
  control = list(factr = 1e3) # Smaller factor means more steps, the default value
  # is 1e7. This is multiplied by .Machine$double.eps
)

# How are we doing on the prediction accuracy ? --------------------------------
joint_mcmc <- readRDS("novel-predictions/mcmc_partially_annotated_no_prior.rds")

pscores_joint <- prediction_score(window(joint_mcmc, start = 5000))
maes_joint    <- 1.0 - sapply(pscores_joint, "[[", "obs") # Getting MAEs

pscores_mle  <- lapply(mle_per_aspect, prediction_score)
maes_mle     <- lapply(pscores_mle, function(i) 1.0 - sapply(i, "[[", "obs"))

pscores_mcmc <- lapply(lapply(mcmc_per_aspect, window, start = 30000), prediction_score)
maes_mcmc    <- lapply(pscores_mcmc, function(i) 1.0 - sapply(i, "[[", "obs"))

# Comparing distributions
data_4_boxplot <- rbind(
  data.frame(Aspect = "Bio. Proc.", Method = "Pooled MCMC", MAE = maes_joint[classes == "biological_process"]),
  data.frame(Aspect = "Bio. Proc.", Method = "Asp. MCMC", MAE = maes_mcmc$biological_process),
  data.frame(Aspect = "Bio. Proc.", Method = "Asp. MLE", MAE = maes_mle$biological_process),
  data.frame(Aspect = "Mol. Fun.", Method = "Pooled MCMC", MAE = maes_joint[classes == "molecular_function"]),
  data.frame(Aspect = "Mol. Fun.", Method = "Asp. MCMC", MAE = maes_mcmc$molecular_function),
  data.frame(Aspect = "Mol. Fun.", Method = "Asp. MLE", MAE = maes_mle$molecular_function),
  data.frame(Aspect = "Cel. Comp.", Method = "Pooled MCMC", MAE = maes_joint[classes == "cellular_component"]),
  data.frame(Aspect = "Cel. Comp.", Method = "Asp. MCMC", MAE = maes_mcmc$cellular_component),
  data.frame(Aspect = "Cel. Comp.", Method = "Asp. MLE", MAE = maes_mle$cellular_component)
)

library(ggplot2)
ggplot(data = data_4_boxplot, aes(y = MAE, x = Aspect)) +
  geom_boxplot(aes(fill = Method)) +
  theme_grey() +
  scale_fill_grey() +
  ylab("Mean Absolute Error (MAE)")

ggsave("novel-predictions/empirical_bayes_mean_square_error.pdf", width = 8, height = 8)

# Comparing coefficients
coef_mle_aspect  <- lapply(mle_per_aspect, coef)
coef_mcmc_aspect <- lapply(lapply(mcmc_per_aspect, window, start = 30000), coef)
coef_joint       <- coef(window(joint_mcmc, start = 5000))

pnames <- names(coef_joint)
data_4_coefplot <- rbind(
  data.frame(Parameter = pnames, Value = coef_joint, Method = "Pooled MCMC", Aspect = "All"),
  data.frame(Parameter = pnames, Value = coef_mle_aspect$biological_process, Method = "Asp. MLE", Aspect = "Bio. Proc."),
  data.frame(Parameter = pnames, Value = coef_mle_aspect$molecular_function, Method = "Asp. MLE", Aspect = "Mol. Fun."),
  data.frame(Parameter = pnames, Value = coef_mle_aspect$cellular_component, Method = "Asp. MLE", Aspect = "Cel. Comp."),
  data.frame(Parameter = pnames, Value = coef_mcmc_aspect$biological_process, Method = "Asp. MCMC", Aspect = "Bio. Proc."),
  data.frame(Parameter = pnames, Value = coef_mcmc_aspect$molecular_function, Method = "Asp. MCMC", Aspect = "Mol. Fun."),
  data.frame(Parameter = pnames, Value = coef_mcmc_aspect$cellular_component, Method = "Asp. MCMC", Aspect = "Cel. Comp.")
)

# Creating the plot
parameter_labels <- c(
  mu_d0  = "mu[paste(d0,1)]",
  mu_d1  = "mu[d10]",
  mu_s0  = "mu[paste(s0,1)]",
  mu_s1  = "mu[s10]",
  psi0 = "psi[paste(0,1)]",
  psi1 = "psi[10]",
  # eta0 = "eta[0]",
  # eta1 = "eta[0]",
  Pi   = "pi"
)

data_4_coefplot$Parameter <- parameter_labels[data_4_coefplot$Parameter]

ggplot(data = data_4_coefplot, aes(y = Value, x=Aspect:Method)) +
  geom_col(aes(fill = Aspect), col = "black", position="dodge") +
  facet_wrap( ~ Parameter, labeller = labeller(Parameter = label_parsed)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  theme(text = element_text(size = 15)) +
  ylab("Parameter Estimate") +
  labs(
    title = "Differences in parameter estimates",
    subtitle = "Comparing Method (Joint or per aspect + MCMC vs MLE)"
    )
  
ggsave("novel-predictions/empirical_bayes_coefs.pdf", width = 8, height = 8)


stop("endofit")

map_estimates <- do.call(rbind, lapply(map_estimates, coef))
ab <- vector("list", ncol(map_estimates))
for (i in 1:ncol(map_estimates)) {
  ab[[i]] <- MASS::fitdistr(
    map_estimates[, i], densfun = "beta",
    start = list(shape1 = ALPHAS[i], shape2 = BETAS[i]),
    lower = .001
  )
}



op <- par(mfrow = c(3, 3))
for (i in 1:length(ab))  {
  curve(dbeta(x, ab[[i]]$estimate[1], ab[[i]]$estimate[2]),
        main = names(ALPHAS)[i], sub = sprintf(
          "MLE: %.4f, Prior: %.4f",
          with(ab[[i]], estimate[1]/sum(estimate)),
          ALPHAS[i]/(ALPHAS[i] + BETAS[i])
        ),
        xlab = "Probability",
        ylab = "Density"
        )
  
  curve(dbeta(x, ALPHAS[i], BETAS[i]),
        main = names(ALPHAS)[i], add = TRUE, col = "red")
}
plot.new()
plot.window(c(0,1),c(0,1))
legend(
  "center",
  cex = 1,
  col = c("black", "red"),
  legend = c(
    "Density based on MLE\nestimates of alpha and beta",
    "Density based on priors\nof alpha and beta"
    )
)

par(op)


# Hierarchical MCMC ------------------------------------------------------------

# Constraint matrix
constr <- kronecker(diag(Ntrees(trees) + 2), matrix(TRUE, nrow = 7, ncol = 7))
constr[nrow(constr) - 1:14 + 1,] <- TRUE
constr <- lower.tri(constr)*constr
diag(constr) <- TRUE
Matrix::image(as(constr, "dgCMatrix"))
Matrix::image(as(crossprod(t(constr)), "dgCMatrix"))

# We want the hyper parameters to move more freely during the warmup stage
Sigma <- diag(c(rep(.001, 7 * 3), rep(0.5, 7 * 2)))

k_ram_wo_sigma <- fmcmc::kernel_ram(
  lb     = .001,
  warmup = 1e3,
  eps    = .001,
  # Sigma  = Sigma,
  ub     = c(rep(.999, 7 * 3), rep(1000, 7 * 2)) #,
  # fixed  = c(rep(FALSE, 7*Ntrees(trees)), rep(TRUE, 7*2)) #,
  # constr = constr
)

k_ram_w_sigma <- fmcmc::kernel_ram(
  lb     = .001,
  warmup = 1e3,
  eps    = .001,
  Sigma  = Sigma,
  ub     = c(rep(.999, 7 * 3), rep(1000, 7 * 2)) #,
  # fixed  = c(rep(FALSE, 7*Ntrees(trees)), rep(TRUE, 7*2)) #,
  # constr = constr
)

k_am <- fmcmc::kernel_am(
  lb     = .001,
  warmup = 5e3,
  eps    = .001, 
  # Sigma  = Sigma,
  ub     = c(rep(.999, 7 * 3), rep(1000, 7 * 2))
)

# Estimating empirical bayes using map
ALPHAS <- c(psi0 = 2, psi1 = 2, mu_d0 = 18, mu_d1 = 10, mu_s0 = 2, mu_s1 = 2, Pi = 2)
BETAS  <- c(psi0 = 18, psi1 = 18, mu_d0 = 2, mu_d1 = 10, mu_s0 = 18, mu_s1 = 18, Pi = 18)

set.seed(123)
# debug(k_ram$proposal)
# debug(aphylo_hier)
# debug()
ans_wo_sigma <- aphylo_hier(
  trees ~ psi + mu_d + mu_s + Pi,
  classes      = as.factor(classes),
  params       = ALPHAS/(ALPHAS + BETAS),
  nsteps       = 2e4,
  kernel       = k_ram_wo_sigma,
  thin         = 1,
  nchains      = 4,
  use_optim    = FALSE,
  # conv_checker = fmcmc::convergence_gelman(freq = 1e3, threshold = 1.001) #,
  hyper_params = cbind(alpha = ALPHAS, beta = BETAS)
)

saveRDS(ans_wo_sigma, "novel-predictions/empirical_bayes_without_sigma.rds")

rm(ans_wo_sigma)
gc(full = TRUE)

set.seed(123)
# debug(k_ram$proposal)
# debug(aphylo_hier)
# debug()
ans_w_sigma <- aphylo_hier(
  trees ~ psi + mu_d + mu_s + Pi,
  classes      = as.factor(classes),
  params       = ALPHAS/(ALPHAS + BETAS),
  nsteps       = 2e4,
  kernel       = k_ram_w_sigma,
  thin         = 1,
  nchains      = 4,
  use_optim    = FALSE,
  # conv_checker = fmcmc::convergence_gelman(freq = 1e3, threshold = 1.001) #,
  hyper_params = cbind(alpha = ALPHAS, beta = BETAS)
)

saveRDS(ans_w_sigma, "novel-predictions/empirical_bayes_with_sigma.rds")
