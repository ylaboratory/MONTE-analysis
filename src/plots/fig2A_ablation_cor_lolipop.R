rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)

compared_metric <- "CPE"
df <- read.csv(paste0("data/monte_outputs/pancancer/monte_ablation_correlation_results.csv"))
df_test_meta <- read.csv("data/metadata/test_pan-cancer_meta.csv")

cancer_num_df <- df_test_meta %>%
  count(df_test_meta$Cancer.type, sort = TRUE)

df$Metric <- NULL
df$Cor_OLS_TopN <- NULL
df$Cor_Var_Weight_TopN <- NULL
colnames(df) <- c("cancer", "OLS", "EB-moderated OLS", "EB + SNR-weighted", "MONTE (EB + SNR + top-N)")

df_long <- df %>%
  pivot_longer(
    cols = -cancer,
    names_to = "method",
    values_to = "value"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("OLS",
                 "EB-moderated OLS",
                 "EB + SNR-weighted",
                 "MONTE (EB + SNR + top-N)")
    )
  )

df_long$cancer <- factor(df_long$cancer, levels=cancer_num_df$`df_test_meta$Cancer.type`)

pdf("results/fig2A_monte_ablation.pdf",width=15, height=5)
ggplot(df_long,
       aes(x = cancer, y = value, color = method)) +
  
  geom_linerange(
    aes(ymin = 0, ymax = value),
    position = position_dodge(width = 0.6),
    linewidth = 0.9
  ) +
  
  geom_point(
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  
  scale_color_manual(
    values = c(
      "OLS" = "#1b9e77",
      "EB-moderated OLS" = "#d95f02",
      "EB + SNR-weighted"="#7570b3",
      "MONTE (EB + SNR + top-N)" = "#296374"
    )
  ) +
  
  labs(
    x = "Cancer type",
    y = paste0("Pearson correlation\n(Predicted purity vs. ", compared_metric, ")"),
    color = "Method"
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size=14, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    plot.subtitle = element_text(color="#57595B"),
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    legend.position = "bottom"
  )

dev.off()
