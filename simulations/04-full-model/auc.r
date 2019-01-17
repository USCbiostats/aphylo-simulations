
library(aphylo)
library(AUC)

source("simulations/00-auc.r")

# dgp  <- readRDS("simulations/dgp.rds")
pred_right <- readRDS("simulations/04-full-model/mcmc_right_prior_prediction.rds")
pred_wrong <- readRDS("simulations/04-full-model/mcmc_wrong_prior_prediction.rds")

ans_right <- aphylo_auc(pred_right)
ans_wrong <- aphylo_auc(pred_wrong)

ans <- rbind(
  data.frame(Prior = "Right", auc = ans_right),
  data.frame(Prior = "Wrong", auc = ans_wrong)
)


# Calculating median 

library(ggplot2)
library(magrittr)
ans %>%
  ggplot(aes(y=auc, x=Prior)) +
  geom_boxplot() +
  geom_hline(yintercept=.95, linetype=2) +
  coord_cartesian(clip="off") +
  annotate("text", x=.5, y=.975, label="0.95", just="left") +
  labs(title = "Distribution of AUC values") + 
  ggsave("simulations/04-full-model/auc.pdf", device = "pdf", width=6, height=6)
