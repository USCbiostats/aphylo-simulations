#!/bin/sh
#SBATCH --account=lc_pdt
#SBATCH --partition=thomas
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB

rm(list = ls())

source("00-gelman-functions.r")

# Checkingout Right prior estimates --------------------------------------------
ans_right_prior <- readRDS("01-gold-standard/mcmc_right_prior.rds")
ans_wrong_prior <- readRDS("01-gold-standard/mcmc_wrong_prior.rds")

ans_right_prior <- ans_right_prior[!sapply(ans_right_prior, inherits, what = "error")]
ans_wrong_prior <- ans_wrong_prior[!sapply(ans_wrong_prior, inherits, what = "error")]

# Extracting the data
N <- length(ans_right_prior)
gelmans_right <- lapply(1:N, get_gelman, obj = ans_right_prior)
gelmans_wrong <- lapply(1:N, get_gelman, obj = ans_wrong_prior)

head(gelmans_right)
head(gelmans_wrong)

gelmans_right <- do.call(rbind, gelmans_right)
gelmans_wrong <- do.call(rbind, gelmans_wrong)

library(ggplot2)
library(magrittr)

gelmans_right <- as.data.frame(gelmans_right)
gelmans_wrong <- as.data.frame(gelmans_wrong)

gelmans_right$Prior <- "Right"
gelmans_wrong$Prior <- "Wrong"

gelmans <- rbind(gelmans_right, gelmans_wrong)

# Saving the data
saveRDS(gelmans, file = "01-gold-standard/gelman.rds")

library(tidyr)

graphics.off()
pdf("01-gold-standard/gelman.pdf")
gelmans %>%
  gather("Statistic", "Score", -Prior) %>%
  subset(Statistic == "mpsrf") %>%
  ggplot(aes(x=Prior, y=Score)) +
  geom_jitter(aes(Prior, Score), alpha=.5, colour="tomato") +
  scale_y_log10() +
  labs(
    title    = "Distribution of Multivariate Gelman Diagnostic",
    subtitle = sprintf(
      "Only %0.2f of the chains did not converged",
      mean(gelmans$mpsrf > 1.1, na.rm = TRUE)
    )
  )
dev.off()

