rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)

df <- read_csv("data/monte_outputs/pancancer/monte_pancancer_model_probe_subset.csv")

df <- df %>%
  mutate(subset_size = as.numeric(subset_size))

df_sum <- df %>%
  group_by(subset_size) %>%
  summarise(
    mse_mean  = mean(test_MSE, na.rm = TRUE),
    mse_sd    = sd(test_MSE, na.rm = TRUE),
    corr_mean = mean(test_correlation, na.rm = TRUE),
    corr_sd   = sd(test_correlation, na.rm = TRUE),
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
    x = "Probe subset size",
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y.left = element_text(size = 16, color = "#E69F00", margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 16, color = "#0072B2", margin = margin(l = 10)),
    axis.text = element_text(size = 12),
    axis.text.y.left = element_text(color = "#E69F00"),
    axis.text.y.right = element_text(color = "#0072B2"),
    panel.grid.minor = element_blank()
  ) + 
  scale_x_log10(
    labels = function(x) {
      parse(text = paste0("10^", log10(x)))
    }
  ) +

pdf("results/fig2C_monte_probe_subset.pdf", width=5, height=4)
print(g)
dev.off()
