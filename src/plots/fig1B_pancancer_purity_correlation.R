rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)
library(patchwork)
library(rlang)
df <- read.csv("data/monte_outputs/pancancer/pancancer_meta_test_with_predictions.csv")

plot_purity <- function(df, metric) {
  metric_sym <- sym(metric)
  
  df_plot <- df %>%
    filter(!is.na(!!metric_sym), !is.na(pred_ESTIMATE))
  
  cor_val <- cor(
    df_plot %>% pull(!!metric_sym),
    df_plot$pred_ESTIMATE,
    method = "pearson"
  )
  
  ggplot(df_plot, aes(x = !!metric_sym, y = pred_ESTIMATE)) +
    geom_point(color = "#547792") +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed",
      color = "#25343F",
      linewidth = 0.7
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0("r = ", round(cor_val, 3)),
      hjust = -0.2, vjust = 2, size = 5
    ) +
    labs(
      x = paste0("True Purity (", metric, ")"),
      y = "Predicted Purity"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text = element_text(size = 12)
    ) +
    coord_equal()
}

g_estimate <- plot_purity(df, "ESTIMATE")
g_cpe <- plot_purity(df, "CPE")
g_lump <- plot_purity(df, "LUMP")

g <- g_estimate + g_lump + g_cpe
print(g)

pdf("results/fig1B_pancancer_purity_correlation.pdf", width=12, height=4)
print(g)
dev.off()
