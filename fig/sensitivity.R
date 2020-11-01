library(aphylo)

estimates <- readRDS("parameter-estimates/mcmc_partially_annotated_no_prior.rds")

# Applying burnin
estimates <- window(
  estimates,
  floor((end(estimates$hist) - start(estimates$hist))/2)
)

# We will iterate through different sets of parameters to investigate what is
# the marginal effect of changing a parameter
N    <- 5
pseq <- seq(.0001, .9999, length.out = N)

# Empty
signif <- vector("list", N)
aucs   <- signif
pscore <- signif

# General output
ans <- vector("list", length(coef(estimates)))
names(ans) <- names(coef(estimates))

for (v in names(ans)) {

  # Capturing current state of the estimates
  params. <- coef(estimates)
  ans[[v]]$signif <- signif
  ans[[v]]$aucs   <- aucs
  ans[[v]]$pscore <- pscore
  
  for (i in 1:N) {
    params.[v] <- pseq[i]
    ps. <- prediction_score(
      estimates,
      loo    = TRUE,
      params = params., 
      ncores = 3L
      ) 
    
    ans[[v]]$signif[[i]] <- sapply(ps., "[[", "pval")
    ans[[v]]$aucs[[i]]   <- sapply(ps., function(i) i$auc$auc)
    ans[[v]]$pscore[[i]] <- 1 - sapply(ps., "[[", "obs")
    
    message("Variable ", v, " test ", i, " done.")
    
  }
    
}

# How much we loose from droping data ------------------------------------------
lose_seq <- seq(0, .9, length.out = 5)
set.seed(123)
ans_n_ann <- vector("list", length(lose_seq))

# Which to estimate
ids <- lapply(estimates$dat, function(d) which(d$tip.annotation[,1] != 9))
for (i in seq_along(lose_seq)) {
  
  estimates_tmp <- estimates
  
  # Available proportion
  estimates_tmp$dat <- lapply(
    estimates_tmp$dat, function(d) {
      
      rel <- which(d$tip.annotation[,1] != 9)
      rel <- sample(rel, length(rel)*lose_seq[i])
      d$tip.annotation[rel, 1] <- 9L
      d
      
    }
    )
  
  estimates_tmp$dat <- do.call(c, estimates_tmp$dat)
  
  if (i == 1) {
    ans_n_ann[[i]] <- prediction_score(
      estimates,
      loo    = TRUE,
      ncores = 3L
    )
  } else {
    
    # Getting the corresponding estimates
    ans_n_ann[[i]] <- predict_pre_order(
      estimates_tmp,
      ids = ids
    )
    
    ans_n_ann[[i]] <- prediction_score(
      x        = estimates_tmp, # ans_n_ann[[i]][[i]][ids[[i]],,drop=FALSE],
      expected = lapply(ans_n_ann[[1]], "[[", "expected") #[ids[[i]],,drop=FALSE]
      )
    
  }
  
  message("Missing ", lose_seq[i], " done.")
  
}


# Is it setting psi to zero help? ----------------------------------------------
params. <- coef(estimates)

ps_baseline <- prediction_score(
  estimates,
  loo = TRUE,
  params = params.,
  ncores = 3L
)

params.[1:2] <- 0

ps_psi_to_0 <- prediction_score(
  estimates,
  loo = TRUE,
  params = params.,
  ncores = 3L
)

auc_baseline <- sapply(ps_baseline, function(i) i$auc$auc)
auc_psi_to_0 <- sapply(ps_psi_to_0, function(i) i$auc$auc)

t.test(
  auc_baseline, 
  auc_psi_to_0,
  alternative = "two.sided",
  paired = TRUE
)

# How about data missing ?
mae_baseline  <- 1- sapply(ps_baseline, "[[", "obs")
mae_drop_data <- lapply(ans_n_ann, function(a) 1 - sapply(a, "[[", "obs"))[-1]
data_missing <- lapply(mae_drop_data, function(m) {
  t.test(m, mae_baseline, paired = TRUE, alternative = "two.sided")
})

data_missing <- do.call(rbind, lapply(data_missing, function(d) {
  data.frame(
    `95% C.I.` = sprintf("[%.2f, %.2f]", d$conf.int[1], d$conf.int[2]),
    `t-stat`   = d$statistic,
    `p-value`  = d$p.value
    )
}))
data_missing$`Prop. drop` <- lose_seq[-1]

print(
  xtable::xtable(data_missing),
  include.rownames = FALSE,
  booktabs = TRUE
)

# For now, saving all
save.image("fig/sensitivity.rda", compress = FALSE)

pseq_nam <- sprintf("%.2f", pseq)

expressions <- list(
  expression(psi[paste(0,1)]),
  expression(psi[10]),
  expression(mu[d01]),
  expression(mu[d10]),
  expression(mu[s01]),
  expression(mu[s10]),
  expression(pi)
)

ylims <- c(.05,.85)
abcut <- .5

graphics.off()
pdf("fig/sensitivity.pdf", width = 6, height = 8)
op <- par(mfrow = c(4, 2), mar = c(2.5,2,2.5,1), oma = c(4,2,0,0))
alldat <- NULL
for (i in seq_along(ans)) {
  # Coloring the closest one
  closest <- which.min(abs(coef(estimates)[i] - pseq))
  colors. <- rep("gray", length(pseq))
  colors.[closest] <- "white"
  
  dat <- do.call(cbind, ans[[i]]$pscore) # - (1 - sapply(ps_baseline, "[[", "obs"))
  
  # Rep sequence
  rep_seq <- rep(nrow(dat), length(pseq_nam))
  alldat <- rbind(
    alldat,
    data.frame(
      MAE    = as.vector(dat),
      par    = names(coef(estimates))[i],#gsub("^expression\\(|\\)$", "", deparse(expressions[[i]])),
      value  = rep(pseq_nam, rep_seq),
      colors = rep(colors., rep_seq)
    )
  )
  
  boxplot(
    x       = dat,
    main    = expressions[[i]],
    names   = pseq_nam, 
    border  = gray(.4),
    col     = colors.,
    outline = FALSE,
    boxwex  = .7,
    ylim    = ylims
  )

  abline(h = abcut, lwd = 1, lty="dashed")
  
}

# Adding plot on the effect of missing data
dat <- do.call(
  cbind,
  lapply(ans_n_ann, function(i) sapply(i, function(j) 1 - j$obs) )
) # - (1 - sapply(ps_baseline, "[[", "obs"))

# par(mar = par()$mar * c(1, 3, 1, 1))

boxplot(
  dat,
  main   = "MAE vs Missing data",
  names  = lose_seq, 
  border = gray(.4),
  outline = FALSE,
  boxwex = .7, 
  ylim   = ylims
)

abline(h = abcut, lwd = 1, lty="dashed")  

par(op)
title(ylab = "Mean Absolute Error (MAE)", xlab = "Value of the parameter")
dev.off()


library(ggplot2)
library(ggridges)

# Creating the plot
parameter_labels <- c(
  psi0   = "psi[paste(0, 1)]",
  psi1   = "psi[10]",
  mu_d0  = "mu[d01]",
  mu_d1  = "mu[d10]",
  mu_s0  = "mu[s01]",
  mu_s1  = "mu[s10]",
  # eta0   = "eta[0]",
  # eta1   = "eta[0]",
  Pi     = "pi"
)

alldat$parlabs <- parameter_labels[alldat$par]

alldat$parlabs <- factor(alldat$parlabs, parameter_labels)

ggplot(data = alldat, mapping = aes(x=MAE, y=value)) +
  geom_density_ridges(aes(fill = colors)) +
  coord_flip() +
  facet_wrap(~parlabs, ncol = 2, labeller = labeller(parlabs = label_parsed)) +
  theme_minimal(base_family = "serif") +
  geom_vline(xintercept = 0.5, lty=2) 

theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  # Adding an horizontal line, and spliting my % of missings
  facet_grid(Prior ~ parameter, labeller = labeller(parameter=label_parsed)) +
  coord_flip()

