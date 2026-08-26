rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)
library(patchwork)

df <- read.csv("data/monte_outputs/pancancer/monte_pancancer_meta_with_predictions_test.csv")
LUMP <- "LUMP"

plot_purity <- function(df, metric, metric2) {
  library(dplyr)
  library(ggplot2)
  
  df_plot <- df %>%
    filter(!is.na(.data[[metric]]), !is.na(.data[[metric2]]))
  
  cor_val <- cor(
    df_plot[[metric]],
    df_plot[[metric2]],
    method = "pearson"
  )
  
  mse_val <- mean(
    (df_plot[[metric]] - df_plot[[metric2]])^2,
    na.rm = TRUE
  )
  print(mse_val)
  
  ggplot(df_plot, aes(x = .data[[metric2]], y = .data[[metric]])) +
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
      x = paste0("Reference purity (", metric2, ")")
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

g_absolute <- plot_purity(df, "ABSOLUTE", "LUMP")
g_estimate <- plot_purity(df, "ESTIMATE", "LUMP")


g <- g_absolute + g_estimate +
  patchwork::plot_layout(ncol = 2) &
  theme(plot.margin = margin(5, 10, 5, 10))

pdf("results/supp_reference_metric_comparison.pdf", width=8, height=4)
print(g)
dev.off()
