rm(list=ls())

library(ggplot2)
library(ggthemes)
library(patchwork)

color_maps <- c(
  "MONTE"="#0A97B0","PAMES"="#F9B487","PureBeta"="#A3B087","InfiniumPurify"="#E8DFCA"
)

# PAMES, PureBeta, InfiniumPurify results
pretrain_runtime <- read.csv("data/benchmark/run_times_pretrained_per-cancer.csv", row.names = 1)
fitted_runtime <- read.csv("data/benchmark/run_times_dataset-fitted_total.csv", row.names = 1)

# MONTE results
monte_training_runtime <- read.csv("data/monte_outputs/pancancer/monte_pancancer_training_time.csv")
monte_training_runtime$n_samples <- NULL
monte_training_runtime$Training.Time..seconds. <- monte_training_runtime$Training.Time..seconds. / 60
colnames(monte_training_runtime) <- colnames(fitted_runtime)
fitted_runtime <- rbind(fitted_runtime, monte_training_runtime)


monte_predict_runtime <- read.csv("data/monte_outputs/pancancer/monte_pancancer_prediction_time_by_cancer.csv")
monte_predict_runtime$prediction_time <- monte_predict_runtime$prediction_time / 60
colnames(monte_predict_runtime) <- c("Cancer", "minutes", "n_samples")
monte_predict_runtime$Method <- "MONTE"
pretrain_runtime <- rbind(pretrain_runtime, monte_predict_runtime)

# merge the results
pretrain_runtime$Method <- factor(pretrain_runtime$Method, levels=c("MONTE", "PAMES", "PureBeta", "InfiniumPurify"))
g_1 <- ggplot(pretrain_runtime) + 
  labs(x="Test samples (per cancer)", y=expression("Runtime (min, " * log[10] * " scale)"),
       title="Runtime for purity estimation \n for pretrained models") +
  geom_point(aes(x=n_samples,y=minutes,color=Method), size=2.2, alpha=0.8) + 
  geom_line(aes(x=n_samples,y=minutes,color=Method), size=1.2, alpha=0.8) +
  scale_y_log10() +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust=0.5),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size=12),
    legend.title = element_text(size=14)
  ) +
  guides(color = guide_legend(nrow = 2)) +
  scale_color_manual(values=color_maps)

pdf("results/fig3B_runtime_predicting.pdf", width=5.5, height=5)
g_1
dev.off()

 