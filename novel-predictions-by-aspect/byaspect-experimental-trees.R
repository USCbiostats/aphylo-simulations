library(aphylo)

trees <- readRDS("data/candidate_trees.rds")

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

go_terms_info <- data.table::fread("data/go_terms_info.csv")

classes <- match(go_terms, go_terms_info$id)
classes <- go_terms_info$aspect[classes]
classes[122] <- "molecular_function" 

if (anyNA(classes))
  stop("One or more GO-terms has no aspect. Need to download the data again!")


# Subset of annotations for biological process
write.csv(data.frame(
  go=sapply(trees[classes == "biological_process"], function(i) colnames(i$tip.annotation)),
  tree=tnames[classes == "biological_process"],
  row.names = NULL
), file = "novel-predictions-by-aspect/bio-process.csv", row.names = FALSE, quote = FALSE)


# Estimating empirical bayes using map
ALPHAS <- c(psi0 = 2, psi1 = 2, mu_d0 = 18, mu_d1 = 10, mu_s0 = 2, mu_s1 = 2, Pi = 2)
BETAS  <- c(psi0 = 18, psi1 = 18, mu_d0 = 2, mu_d1 = 10, mu_s0 = 18, mu_s1 = 18, Pi = 18)

bpriors <- bprior(shape1 = ALPHAS, shape2 = BETAS)

# How evenly distributed are the types?
prop_types <- sapply(trees, function(i) mean(with(i, c(tip.type, node.type)[reduced_pseq])))
tapply(X = prop_types, INDEX = classes, mean)

# How about the number of 0s and 1s
zeros_and_ones <- sapply(trees, function(i) aphylo:::fast_table_using_labels(i$tip.annotation, c(0, 1, 9)))
zeros_and_ones <- t(zeros_and_ones)
zeros_and_ones <- tapply(
  X = asplit(zeros_and_ones, 1),
  INDEX = classes,
  function(i) {
    data.frame(
      Avg    = colMeans(do.call(rbind, i)),
      Counts = colSums(do.call(rbind, i))
    )
    }
  )

zeros_and_ones <- lapply(names(zeros_and_ones), function(i) {
  cbind(Aspect = i, Annotation = c(0, 1, NA), zeros_and_ones[[i]])
})
zeros_and_ones <- do.call(rbind, zeros_and_ones)
print(zeros_and_ones, digits = 3)


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

saveRDS(
  mle_per_aspect,
  "novel-predictions-by-aspect/byaspect-experimental-trees_mle.rds"
  )

saveRDS(
  mcmc_per_aspect,
  "novel-predictions-by-aspect/byaspect-experimental-trees_mcmc.rds"
  )


# Biological Process the MCMC is doing way better, so it is a good idea to
# rerun the MLE starting from the MCMC estimates
mle_per_aspect$biological_process <- aphylo_mle(
  trees[classes == "biological_process"] ~ psi + mu_d + mu_s + Pi,
  params = coef(window(mcmc_per_aspect$biological_process, start = 30000)),
  
  # Smaller factor means more steps, the default value
  # is 1e7. This is multiplied by .Machine$double.eps
  control = list(factr = 1e3) 
  
)

# How are we doing on the prediction accuracy ? --------------------------------
joint_mcmc <- readRDS("novel-predictions/mcmc_partially_annotated_no_prior.rds")

pscores_joint <- prediction_score(window(joint_mcmc, start = 5000))
maes_joint    <- 1.0 - sapply(pscores_joint, "[[", "obs") # Getting MAEs

pscores_mle  <- lapply(mle_per_aspect, prediction_score)
maes_mle     <- lapply(pscores_mle, function(i) 1.0 - sapply(i, "[[", "obs"))

pscores_mcmc <- lapply(lapply(mcmc_per_aspect, window, start = 30000L), prediction_score)
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

ggsave(
  "novel-predictions-by-aspect/byaspect-experimental-trees_mae.pdf",
  width = 8, height = 8
  )

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
  
ggsave("novel-predictions-by-aspect/byaspect-experimental-trees_coefs.pdf", width = 8, height = 8)

# Tabulating the results -------------------------------------------------------

# Accuracy scores: For this we take each set of predictions as a whole,
# instead of by tree, I mean, unlisting all leafs

# For the by aspect
maes_mcmc <- lapply(pscores_mcmc, function(ps) {
  dat. <- do.call(
    rbind,
    lapply(ps, function(p) cbind(p$predicted, p$expected))
    )
  
  idx <- which(dat.[,2] != 9L)
  
  cbind(
    AUC = auc(pred = dat.[,1], labels = dat.[,2])$auc,
    MAE = mean(abs(dat.[idx,2] - dat.[idx,1]))
  )
})

maes_mcmc <- t(do.call(rbind, maes_mcmc))
colnames(maes_mcmc) <- names(pscores_mcmc)
maes_mcmc[] <- sprintf("%.2f", maes_mcmc)

# for the joint case (need to match the names first)
names_in_joint <- sapply(pscores_joint, function(i) colnames(i$predicted))

classes_joint <- match(names_in_joint, go_terms_info$id)
classes_joint <- go_terms_info$aspect[classes_joint]
classes_joint[122] <- "molecular_function" # For some reason is NA

# Re-grouping by aspect
maes_joint <- unique(classes_joint)
names(maes_joint) <- maes_joint
maes_joint <- lapply(maes_joint, function(c.) {
  pscores_joint[which(classes_joint == c.)]
})

# Computing AUCs and MAEs
maes_joint <- lapply(maes_joint, function(ps) {
  dat. <- do.call(
    rbind,
    lapply(ps, function(p) cbind(p$predicted, p$expected))
  )
  
  idx <- which(dat.[,2] != 9L)
  
  cbind(
    AUC = auc(pred = dat.[,1], labels = dat.[,2])$auc,
    MAE = mean(abs(dat.[idx,2] - dat.[idx,1]))
  )
})

maes_joint <- t(do.call(rbind, maes_joint))
colnames(maes_joint) <- colnames(maes_mcmc)
maes_joint[] <- sprintf("%.2f", maes_joint)

# A table can be neat as well
neat_table <- cbind(
  Pooled = coef(window(joint_mcmc, start = 5000L)),
  sapply(mcmc_per_aspect, function(m) coef(window(m, start = 25000L)))
)

neat_table[] <- sprintf("%.2f", neat_table[])

# Tree count
neat_table <- rbind(
  neat_table,
  Trees = sprintf(
    "%d",
    c(Ntrees(joint_mcmc), sapply(mcmc_per_aspect, Ntrees))
  ),
  cbind(Pooled = "-", maes_mcmc),
  cbind(Pooled = "-", maes_joint)
)

print(xtable::xtable(neat_table), booktabs = TRUE)



# hierarchical model for molecular function ------------------------------------

# Obtaining the alpha and beta
ab <- mcmc_per_aspect$molecular_function$hist
ab <- window(ab, start = 25000)
ab <- do.call(rbind, as.list(ab))

fitted_dist <- vector("list", ncol(ab))
for (i in 1:ncol(ab)) {
  fitted_dist[[i]] <- MASS::fitdistr(
    x = ab[,i],
    densfun = "beta", 
    start   = list(shape1 = 10, shape2 = 10)
    )
  message(i, " done...")
}

alpha_beta <- t(sapply(fitted_dist, coef))
dimnames(alpha_beta) <- list(
  names(ALPHAS),
  c("alpha", "beta")
)

# We want the hyper parameters to move more freely during the warmup stage
n_mol_fun_trees <- sum(classes == "molecular_function")
Sigma <- diag(c(rep(.001, 7 * n_mol_fun_trees), rep(0.5, 7 * 2)))


k_ram <- fmcmc::kernel_ram(
  lb     = .001,
  warmup = 1e3,
  eps    = .001,
  # Sigma  = Sigma,
  ub     = c(rep(.999, 7 * n_mol_fun_trees), rep(1000, 7 * 2)),
  fixed  = c(rep(FALSE, 7 * n_mol_fun_trees), rep(TRUE, 7*2)) #,
  # constr = constr
)

set.seed(123)
# debug(k_ram$proposal)
# debug(aphylo_hier)
# debug()
hierarchical_mol_fun <- aphylo_hier(
  trees[classes == "molecular_function"] ~ psi + mu_d + mu_s + Pi,
  classes      = 1:n_mol_fun_trees,
  params       = coef(window(mcmc_per_aspect$molecular_function, start = 25000)),
  nsteps       = 4e4,
  kernel       = k_ram,
  thin         = 1,
  nchains      = 2,
  use_optim    = FALSE,
  # conv_checker = fmcmc::convergence_gelman(freq = 1e3, threshold = 1.001) #,
  hyper_params = alpha_beta
)

saveRDS(
  hierarchical_mol_fun,
  "novel-predictions-by-aspect/empirical_bayes_molecular_function.rds"
  )

k_ram2 <- fmcmc::kernel_ram(
  lb     = .001,
  warmup = 1e3,
  eps    = .001,
  # Sigma  = Sigma,
  ub     = c(rep(.999, 7 * n_mol_fun_trees), rep(1000, 7 * 2)),
  fixed  = c(rep(FALSE, 7 * n_mol_fun_trees), rep(TRUE, 7*2)) #,
  # constr = constr
)
hierarchical_mol_fun <- readRDS(
  "novel-predictions-by-aspect/empirical_bayes_molecular_function.rds"
  )

set.seed(71623)
hierarchical_mol_fun <- aphylo_hier(
  trees[classes == "molecular_function"] ~ psi + mu_d + mu_s + Pi,
  classes      = 1:n_mol_fun_trees,
  params       = coef(window(mcmc_per_aspect$molecular_function, start = 25000)),
  params0      = colMeans(do.call("rbind", window(hierarchical_mol_fun$hist, start=30000))),
  nsteps       = 4e4,
  kernel       = k_ram2,
  thin         = 1,
  nchains      = 2,
  use_optim    = FALSE,
  conv_checker = fmcmc::convergence_gelman(freq = 500, threshold = 1.001),
  hyper_params = alpha_beta
)

# Looking at prediction accuracy with this

# Halfing
smaller_sample <- window(hierarchical_mol_fun$hist, start = 2e4)
1 - coda::rejectionRate(smaller_sample) # ~ 34% acceptance

# Freeing some space
hierarchical_mol_fun$hist <- NULL

# we start by retrieving parameters
n_mol_fun <- sum(classes == "molecular_function")
individual_parameters <- lapply(1:n_mol_fun, function(i) {
  
  # library(coda)
  
  # Forming the names
  pnames <- sprintf("%s_class%03d", names(ALPHAS), i)
  
  structure(
    summary(smaller_sample[,pnames])$statistics[,"Mean"],
    names = names(ALPHAS)
  )
  
})

individual_gelman <- lapply(1:n_mol_fun, function(i) {
  
  # Forming the names
  pnames <- sprintf("%s_class%03d", names(ALPHAS), i)
  
  coda::gelman.diag(smaller_sample[,pnames])
  
  
})

# Checking that colnames matches
emp_bayes_cnames <- sapply(
  hierarchical_mol_fun$dat,
  function(i) colnames(i$tip.annotation)
)

aspect_cnames <- sapply(
  mcmc_per_aspect$molecular_function$dat,
  function(i) colnames(i$tip.annotation)
)

# These should match
stopifnot(all(emp_bayes_cnames == aspect_cnames))

predictions_emp_bayes <- vector("list", length(emp_bayes_cnames))
for (i in seq_along(predictions_emp_bayes)) {
  
  # Figuring out the ids
  ids <- which(
    mcmc_per_aspect$molecular_function$dat[[i]]$tip.annotation != 9
    )
  
  predictions_emp_bayes[[i]] <- predict(
    mcmc_per_aspect$molecular_function,
    which.tree = i,
    ids        = list(ids),
    params     = individual_parameters[[i]]
  )[[1]][ids,]
  
  predictions_emp_bayes[[i]] <- cbind(
    expected  = mcmc_per_aspect$molecular_function$dat[[i]]$tip.annotation[ids],
    predicted = predictions_emp_bayes[[i]]
  )
    
  message(i, " done.")
}

predictions_emp_bayes <- do.call(rbind, predictions_emp_bayes)

emp_bayes_auc <- auc(
  pred   = predictions_emp_bayes[,2],
  labels = predictions_emp_bayes[,1]
  )

plot(emp_bayes_auc)

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
