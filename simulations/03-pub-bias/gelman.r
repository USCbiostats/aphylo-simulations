rm(list = ls())

source("simulations/00-gelman-functions.r")

# Checkingout Right prior estimates --------------------------------------------
ans_MCMC_right_prior <- readRDS("simulations/03-pub-bias/mcmc_right_prior_estimates.rds")
ans_MCMC_wrong_prior <- readRDS("simulations/03-pub-bias/mcmc_wrong_prior_estimates.rds")

# Extracting the data
N <- length(ans_MCMC_right_prior)
gelmans_right <- do.call(rbind, parallel::mclapply(1:N, get_gelman, obj = ans_MCMC_right_prior, mc.cores=10L))
gelmans_wrong <- do.call(rbind, parallel::mclapply(1:N, get_gelman, obj = ans_MCMC_wrong_prior, mc.cores=10L))

library(ggplot2)
library(magrittr)

gelmans_right <- as.data.frame(gelmans_right)
gelmans_wrong <- as.data.frame(gelmans_wrong)

gelmans_right$Prior <- "Right"
gelmans_wrong$Prior <- "Wrong"

gelmans <- rbind(gelmans_right, gelmans_wrong)

# Saving the data
saveRDS(gelmans, file = "simulations/03-pub-bias/gelman.rds")

library(tidyr)

graphics.off()
pdf("simulations/03-pub-bias/gelman.pdf")
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

