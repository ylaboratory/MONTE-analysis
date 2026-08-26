rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)
library(patchwork)

df_all <- read_csv("data/external_dataset/pcawg_fine_tuning_sample_size_results.csv")
cancer_types <- unique(df_all$cancer)

plot_sample_subset <- function(df, cancer, cancer_name)
{
  df <- df[df$cancer == cancer,]
  df <- df %>%
    mutate(subset_size = as.numeric(train_size))
  
  df_sum <- df %>%
    group_by(subset_size) %>%
    summarise(
      mse_mean  = mean(mse, na.rm = TRUE),
      mse_sd    = sd(mse, na.rm = TRUE),
      corr_mean = mean(pearson_corr, na.rm = TRUE),
      corr_sd   = sd(pearson_corr, na.rm = TRUE),
      .groups = "drop"
    )
  
  scale_factor <- max(df_sum$mse_mean + df_sum$mse_sd, na.rm = TRUE) /
    max(df_sum$corr_mean + df_sum$corr_sd, na.rm = TRUE)
  
  df_sum <- df_sum %>%
    mutate(
      corr_scaled = corr_mean * scale_factor,
      corr_sd_scaled = corr_sd * scale_factor
    )
  
  g <- ggplot(df_sum, aes(x = subset_size)) +
    # MSE ribbon + line
    geom_ribbon(
      aes(
        ymin = mse_mean - mse_sd,
        ymax = mse_mean + mse_sd
      ),
      fill = "#E69F00",
      alpha = 0.18
    ) +
    geom_line(
      aes(y = mse_mean),
      color = "#E69F00",
      linewidth = 1.2
    ) +
    geom_point(
      aes(y = mse_mean),
      color = "#E69F00",
      size = 2
    ) +
    
    # Correlation ribbon + line (scaled)
    geom_ribbon(
      aes(
        ymin = corr_scaled - corr_sd_scaled,
        ymax = corr_scaled + corr_sd_scaled
      ),
      fill = "#0072B2",
      alpha = 0.16
    ) +
    geom_line(
      aes(y = corr_scaled),
      color = "#0072B2",
      linewidth = 1.2
    ) +
    geom_point(
      aes(y = corr_scaled),
      color = "#0072B2",
      size = 2
    ) +
    scale_y_continuous(
      name = "Mean Squared Error",
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Pearson correlation"
      )
    ) +
    labs(
      title = cancer_name,
      x = "Sample subset size",
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      axis.title.y.left = element_text(size = 16, color = "#E69F00", margin = margin(r = 10)),
      axis.title.y.right = element_text(size = 16, color = "#0072B2", margin = margin(l = 10)),
      axis.text = element_text(size = 14),
      axis.text.y.left = element_text(color = "#E69F00"),
      axis.text.y.right = element_text(color = "#0072B2"),
      panel.grid.minor = element_blank()
    )
  return (g)
}

g_ov <- plot_sample_subset(df_all, "ov_au", "PCAWG (OV)")
g_paca <- plot_sample_subset(df_all, "paca_au", "PCAWG (PACA)")
g_prad <- plot_sample_subset(df_all, "prad_ca", "PCAWG (PRAD)")

g <- g_ov + g_paca + g_prad

pdf("results/fig4B_pcawg_sample_subset.pdf", width=15, height=4)
print(g)
dev.off()
