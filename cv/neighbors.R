library(aphylo)

cross_validation <- readRDS("cv/cross_validation.rds")

closests_0s_and_1s <- function(aphylo, expected, predicted) {
  
  # Calculating the distance to each pair
  tokeep <- which(expected != 9L)
  expected  <- expected[tokeep, ,drop=FALSE]
  predicted <- predicted[tokeep, ,drop=FALSE]
  D <- igraph::graph_from_edgelist(aphylo$edge, directed = FALSE)
  D <- igraph::distances(D)
  
  ones  <- rownames(expected)[which(expected == 1L)]
  zeros <- setdiff(rownames(expected), ones)
  ones  <- as.integer(ones)
  zeros <- as.integer(zeros)
  
  # Number of zeros
  diag(D) <- 1e4
  D <- D[as.integer(rownames(expected)), , drop = FALSE]
  
  n <- nrow(D)
  
  # Closests zero
  if (length(zeros)) {
    zero_count <- D[, zeros, drop=FALSE]
    zero_count <- zero_count[cbind(1:n, max.col(-zero_count))]
    zero_count[zero_count == 1e4] <- NA
  } else {
    zero_count <- matrix(NA, nrow=n, ncol=1L)
  }
  
  # Closests one
  if (length(ones)) {
    one_count <- D[, ones, drop=FALSE]
    one_count <- one_count[cbind(1:n, max.col(-one_count))]
    one_count[one_count == 1e4] <- NA
  } else {
    one_count <- matrix(NA, nrow=n, ncol=1L)
  }
  
  data.frame(
    expected      = unname(expected),
    predicted     = unname(predicted),
    edgeid        = rownames(expected),
    term          = colnames(expected),
    closests_zero = zero_count,
    closests_one  = one_count,
    n_zeros       = length(zeros),
    n_ones        = length(ones)
    )

}

current <- 0L
ans <- lapply(cross_validation, function(d) {
  current <<- current + 1L
  cat("Going for ", names(cross_validation)[current], "\n")
  closests_0s_and_1s(
    d$estimates$dat$tree,
    d$pscore$expected,
    d$pscore$predicted
  )
})

ans <- do.call(rbind, ans)
ans$tree <- gsub("[.][0-9]+$", "", rownames(ans))
rownames(ans) <- NULL

readr::write_csv(ans, "cv/neighbors.csv")

# Plots ------------------------------------------------------------------------

library(dplyr)
library(magrittr)
library(ggplot2)

# Closests same vs Closests different ------------------------------------------
same_vs_diff <- 
  ans[complete.cases(ans),] %>%
  as_tibble %>%
    mutate(
      closests_same = if_else(expected == 1, closests_one, closests_zero),
      closests_diff = if_else(expected == 0, closests_one, closests_zero)
    ) %>%
  group_by(closests_same, closests_diff, expected) %>%
  summarise(cell = 1 - mean(abs(expected - predicted))) %>%
  arrange(closests_same, closests_diff)

ggplot(same_vs_diff, aes(x = closests_same, y = closests_diff)) +
  geom_raster(aes(fill = cell)) + 
  facet_grid(
    ~expected,
    labeller = labeller(
      expected = c(`0`="Function: Not present", `1` = "Function: Present")
      )
    ) +
  theme_bw() +
  scale_fill_viridis_c("Level of\nAccuracy") +
  xlab("Distance to the closests\nleaf with the same annotation") +
  ylab("Distance to the closests\nleaf with the the opposite annotation")


image(
  x = ans2$closests_zero,
  y = ans2$closests_one,
  z = ans2$cell
)
