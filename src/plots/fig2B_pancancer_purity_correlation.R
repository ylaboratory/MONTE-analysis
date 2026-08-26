rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)
library(patchwork)

df <- read.csv("data/monte_outputs/pancancer/monte_pancancer_meta_with_predictions_test.csv")
pred_metric <- "predicted_CPE"

plot_purity <- function(df, metric, pred_metric) {
  library(dplyr)
  library(ggplot2)
  
  df_plot <- df %>%
    filter(!is.na(.data[[metric]]), !is.na(.data[[pred_metric]]))
  
  cor_val <- cor(
    df_plot[[metric]],
    df_plot[[pred_metric]],
    method = "pearson"
  )
  
  mse_val <- mean(
    (df_plot[[metric]] - df_plot[[pred_metric]])^2,
    na.rm = TRUE
  )
  print(mse_val)
  
  ggplot(df_plot, aes(x = .data[[pred_metric]], y = .data[[metric]])) +
    geom_point(color = "#547792", size=.8) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed",
      color = "#25343F",
      linewidth = 0.7
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0("rho == ", round(cor_val, 3)),
      parse = TRUE,
      hjust = -0.2, vjust = 2, size = 5
    ) +
    labs(
      y = paste0("Reference purity (", metric, ")"),
      x = paste0("MONTE predicted purity")
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 15, margin = margin(t = 10)),
      axis.title.y = element_text(size = 15, margin = margin(r = 10)),
      axis.text = element_text(size = 12)
    ) +
    coord_equal()
}

g_absolute <- plot_purity(df, "ABSOLUTE", pred_metric)
g_estimate <- plot_purity(df, "ESTIMATE", pred_metric)
g_cpe <- plot_purity(df, "CPE", pred_metric)
# g_lump <- plot_purity(df, "LUMP", pred_metric)

g <- g_cpe + g_absolute + g_estimate +
  patchwork::plot_layout(ncol = 3) &
  theme(plot.margin = margin(5, 10, 5, 10))

pdf("results/fig2B_monte_pancancer_purity_correlation.pdf", width=12, height=4)
print(g)
dev.off()
