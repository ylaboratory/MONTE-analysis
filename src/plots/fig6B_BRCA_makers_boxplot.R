rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)
library(ggpubr)

df <- read.csv("data/application/BRCA_marker_effect_change_summary.csv", check.names = FALSE)
df$`BRCA markers` <- log(df$`BRCA markers`)
df$`Non-markers` <- log(df$`Non-markers`)

df_long <- df %>%
  pivot_longer(
    cols = c(`BRCA markers`, `Non-markers`),
    names_to = "group",
    values_to = "value"
  ) %>%
  drop_na(value)

# wilcox_res <- wilcox.test(
#   df$`BRCA markers`,
#   df$`Non-markers`,
#   paired = TRUE,
#   alternative = "greater"
# )
# 
# p_val <- wilcox_res$p.value
# p_val

pdf("results/fig6B_BRCA_maker_genes.pdf", width=5, height=6.5)
ggplot(df_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = TRUE, alpha = 0.7) +
  # stat_compare_means(
  #   comparisons = list(c("BRCA markers", "Non-markers")),
  #   method = "wilcox.test",
  #   paired = TRUE,
  #   method.args = list(alternative = "greater"),
  #   label = "p.signif",
  #   size=6
  # ) +
  labs(
    x = NULL,
    y = "Effect size change",
    subtitle = "Permutation test p=0.0002"
  ) +
  theme_bw() +
  theme(
    plot.subtitle = element_text(size=14, hjust=0.5),
    legend.position = "none",
    axis.title = element_text(size=16),
    axis.text = element_text(size=14),
  ) +
  scale_fill_manual(values=c("BRCA markers"="#C08552", "Non-markers"="#FFF8F0"))
dev.off()

df_coeff <- read.csv("data/application/BRCA_marker_effect_change_sampling.csv")
df_coeff$log_coefficient_change <- log(df_coeff$coefficient_change)

ggplot(df_coeff, aes(x = log_coefficient_change, fill = type)) +
  geom_density(alpha = 0.5) +
  xlim(0, 2.5) +
  theme_bw() +
  labs(
    x = "Log(Coefficient change)",
    y = "Density",
    fill = "Type"
  )
