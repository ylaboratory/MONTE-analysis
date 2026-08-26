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

inf_file <- 'data/benchmark/infiniumpurify_summary_table.csv'
pure_file <- 'data/benchmark/purebeta_trainval_to_test_summary.csv'
infinium_training_runtime <- read.csv(inf_file)
infinium_training_runtime <- infinium_training_runtime[infinium_training_runtime$method == "refFree",]
infinium_training_runtime <- infinium_training_runtime[!is.na(infinium_training_runtime$cor_vs_CPE), ]
names(infinium_training_runtime)[names(infinium_training_runtime) == "purity_time_sec"] <- "fitted_time_s"

purebeta_training_runtime <- read.csv(pure_file)
purebeta_training_runtime$fitted_time_s <- purebeta_training_runtime$train_time_sec + purebeta_training_runtime$purity_time_sec

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

fitted_runtime$Method <- factor(fitted_runtime$Method, levels=c("MONTE", "PureBeta", "InfiniumPurify"))
# g_2 <- ggplot(fitted_runtime, aes(x = Method, y = minutes, fill = Method)) +
#   geom_col(width=0.7) +
#   labs(
#     x = "",
#     y = "Time (min)",
#     title = "Runtime for purity estimation \n for dataset-fitted models (aggregate)"
#   ) +
#   scale_y_log10() +
#   theme_bw() +
#   theme(
#     plot.title = element_text(size = 16, hjust=0.5),
#     axis.title.x = element_text(size = 16, margin = margin(t = 10)),
#     axis.title.y = element_text(size = 16, margin = margin(r = 10)),
#     axis.text = element_text(size = 14),
#     legend.position = "none"
#   ) +
#   scale_fill_manual(values = color_maps)


method_levels <- c(
  "MONTE",
  "PureBeta",
  "InfiniumPurify"
)

bar_width <- 0.7

monte_runtime <- fitted_runtime %>%
  filter(Method == "MONTE") %>%
  summarise(
    total_minutes = sum(minutes, na.rm = TRUE)
  ) %>%
  mutate(
    Method = "MONTE"
  )

purebeta_project_runtime <- purebeta_training_runtime %>%
  group_by(project) %>%
  summarise(
    minutes = sum(fitted_time_s, na.rm = TRUE) / 60,
    .groups = "drop"
  ) %>%
  filter(
    !is.na(minutes),
    minutes > 0
  ) %>%
  arrange(desc(minutes))


# PureBeta total runtime
purebeta_total <- purebeta_project_runtime %>%
  summarise(
    total_minutes = sum(minutes, na.rm = TRUE)
  ) %>%
  mutate(
    Method = "PureBeta"
  )

infinium_project_runtime <- infinium_training_runtime %>%
  group_by(project) %>%
  summarise(
    minutes = sum(fitted_time_s, na.rm = TRUE) / 60,
    .groups = "drop"
  ) %>%
  filter(
    !is.na(minutes),
    minutes > 0
  ) %>%
  arrange(desc(minutes))


# InfiniumPurify total runtime
infinium_total <- infinium_project_runtime %>%
  summarise(
    total_minutes = sum(minutes, na.rm = TRUE)
  ) %>%
  mutate(
    Method = "InfiniumPurify"
  )


bar_runtime <- bind_rows(
  monte_runtime,
  purebeta_total,
  infinium_total
) %>%
  mutate(
    Method = factor(
      Method,
      levels = method_levels
    ),
    x = match(
      as.character(Method),
      method_levels
    )
  )


purebeta_separators <- purebeta_project_runtime %>%
  mutate(
    boundary = cumsum(minutes),
    total = sum(minutes)
  ) %>%
  
  # Remove the final boundary (= top of the bar)
  filter(boundary < total) %>%
  
  mutate(
    Method = "PureBeta",
    x = match("PureBeta", method_levels),
    xmin = x - bar_width / 2,
    xmax = x + bar_width / 2
  )


infinium_separators <- infinium_project_runtime %>%
  mutate(
    boundary = cumsum(minutes),
    total = sum(minutes)
  ) %>%
  
  # Remove the final boundary (= top of the bar)
  filter(boundary < total) %>%
  
  mutate(
    Method = "InfiniumPurify",
    x = match("InfiniumPurify", method_levels),
    xmin = x - bar_width / 2,
    xmax = x + bar_width / 2
  )


separator_df <- bind_rows(
  purebeta_separators,
  infinium_separators
)

print(bar_runtime)

print(
  purebeta_project_runtime %>%
    summarise(
      sum_project_minutes = sum(minutes)
    )
)

print(
  infinium_project_runtime %>%
    summarise(
      sum_project_minutes = sum(minutes)
    )
)

print(separator_df)


g_2 <- ggplot() +
geom_col(
  data = bar_runtime,
  aes(
    x = x,
    y = total_minutes,
    fill = Method
  ),
  width = bar_width,
  color = "black",
  linewidth = 0.35
) +
geom_segment(
  data = separator_df,
  aes(
    x = xmin,
    xend = xmax,
    y = boundary,
    yend = boundary
  ),
  inherit.aes = FALSE,
  color = "black",
  linewidth = 0.35,
  linetype = "22"
) +
labs(
  x = "",
  y = expression("Runtime (min, " * log[10] * " scale)"),
  title = "Runtime for purity estimation \n for dataset-fitted models (aggregate)"
) +
scale_x_continuous(
  breaks = seq_along(method_levels),
  labels = method_levels
) +
scale_y_log10() +
scale_fill_manual(
  values = color_maps
) +
theme_bw() +
  
  theme(
    plot.title = element_text(
      size = 16,
      hjust = 0.5
    ),
    
    axis.title.x = element_text(
      size = 16,
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      size = 16,
      margin = margin(r = 10)
    ),
    
    axis.text = element_text(
      size = 14
    ),
    
    legend.position = "none"
  )


g_2
pdf("results/fig3C_runtime_training.pdf", width=5.5, height=5)
g_2
dev.off()
