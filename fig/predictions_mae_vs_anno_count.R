library(data.table)

predlll <- fread("parameter-estimates/predictions.csv.gz")

# How well it correlates as a function of # of annotations ---------------------
maes <- predlll[, .(mae = mean(abs(score - truth)), n = .N, n0 = .N - sum(truth)), by = .(term, panther)]
maes[,plot(n0, mae)]

library(ggplot2)

graphics.off()
pdf("fig/predictions_mae_vs_anno_count.pdf", width=6, height = 5)
set.seed(123)
ggplot(maes, aes(x = n, y = mae)) +
  theme_gray() +
  geom_jitter(aes(color = as.factor(n0)), cex = 2) +
  labs(color = "Number of\n\"absent\"\nannotations", y = "MAE", x = "Total Number of\nannotations (log-scale)") +
  scale_x_log10() + 
  geom_smooth() +
  theme_minimal(base_family = "serif") 
dev.off()