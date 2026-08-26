rm(list = ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)
library(patchwork)

plot_pcawg <- function(cancer, cancer_name) {
  df <- read.csv(paste0("data/external_dataset/pcawg_", cancer, "_purity_predictions.csv"))
  colnames(df) <- c("donor_id", "purity", "No", "Yes")
  
  df_long <- df %>% 
    pivot_longer(
      cols = c(Yes, No),
      names_to = "fine_tuned",
      values_to = "predicted_purity"
    )
  
  cor_df <- df_long %>%
    group_by(fine_tuned) %>%
    summarise(
      cor_val = cor(purity, predicted_purity, method = "pearson", use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(fine_tuned, ":~rho==", round(cor_val, 3)),
      x = 0.02,
      y = ifelse(fine_tuned == "Yes", 0.98, 0.88)
    )
  
  cor_df$text_color <- ifelse(
    cor_df$fine_tuned == "No",
    "#57595B",   # darker No text
    "#547792"    # Yes text
  )
  
  ggplot(df_long, aes(x = predicted_purity, y = purity, color = fine_tuned)) + 
    geom_point(alpha=0.8) +
    geom_text(
      data = cor_df,
      aes(x = x, y = y, label = label),
      color = cor_df$text_color,
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 5,
      show.legend = FALSE,
      parse=TRUE
    ) +
    labs(
      y = "PCAWG reference purity",
      x = "MONTE predicted purity",
      title = paste0("PCAWG (", cancer_name,")"),
      color = "With transfer learning"
    ) +
    xlim(0, 1) + ylim(0, 1) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed",
      color = "#57595B",
      linewidth = 0.7
    ) + 
    scale_color_manual(values = c("Yes" = "#547792", "No" = "#BFC6C4")) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 14),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    ) +
    coord_equal()
}

# "ov_au", "paca_au", "prad_ca"
g_ov <- plot_pcawg("ov_au", "OV")
g_paca <- plot_pcawg("paca_au", "PACA")
g_prad <- plot_pcawg("prad_ca", "PRAD")

g <- g_ov + g_paca + g_prad +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf("results/fig4A_pcawg_fine_tune.pdf", width=12, height=5)
print(g)
dev.off()