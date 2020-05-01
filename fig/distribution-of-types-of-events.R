library(aphylo)

estimates <- readRDS("novel-predictions/mcmc_partially_annotated_no_prior.rds")

# Extracting the proportions ---------------------------------------------------
types <- sapply(estimates$dat, function(d) {
  
  # mean(1 - d$node.type)
  aphylo:::fast_table_using_labels(
    d$node.type,
    c(0, 1)
  )
  
})

graphics.off()
pdf("fig/distribution-event-type.pdf", width = 6, height = 4)
op <- par(mar = par("mar") * c(.2, 1,.2, .2), xpd = NA)
barplot(
  types[, order(colSums(types))],
  las    = 3,
  col    = c("black", gray(.4)),
  space  = 0,
  border = "transparent",
  xlab   = "Phylogeny",
  ylab   = "Number of internal nodes",
  main   = NULL
  )
legend(
  x = 1, y = 1000,
  legend = c("Duplication", "Other (mostly speciation)"),
  bty = "n",
  fill = c("black", gray(.5))
)

# Second plot, a histogram
op2 <- par(omd = c(.3, .85, .5, .9), new = FALSE)
plot.window(xlim = c(0,.55), ylim = c(0, 35), new = FALSE)
hist(
  types[1,]/colSums(types),
  add    = TRUE,
  col    = "gray",
  border = gray(.4)
  )

# Adding text
title(
  main = c("Distribution of proportions\nof duplication events"),
  font.main = 1,
  cex.main=.9
  )
title(ylab = "Number of trees", line = 2)

axis(1);axis(2)
par(op);par(op2)

dev.off()

# Extracting the annotation types ----------------------------------------------
labels <- sapply(estimates$dat, function(d) {
  structure(
    aphylo:::fast_table_using_labels(
      d$tip.annotation,
      c(0, 1, 9)
      ),
    names = c("0", "1", "9")
  )
})
labels <- t(labels)[, -3]

labels <- labels[order(rowSums(labels)),]

graphics.off()
pdf("fig/distribution-annotation-type.pdf", width = 6, height = 4)
op <- par(mar = par("mar") * c(.2, 1, .2, .2))
barplot(
  t(labels),
  border = "transparent",
  col    = c("black", gray(.5)),
  space  = 0,
  ylab   = c("Number of annotations", "(log-scale)"),
  xlab   = "Trees",
  las    = 3, axes = FALSE,
  log    = "y"
  )

legend(
  x = 10, y = 40,
  title  = "Annotation type",
  legend = c("Zero", "One"),
  fill   = c("black", gray(.5)),
  bty    = "n"
  )
axis(2, at = c(2, 10, 30, 60))
par(op)
dev.off()