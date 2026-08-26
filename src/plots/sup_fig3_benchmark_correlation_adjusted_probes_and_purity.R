rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)
library(rstatix)
library(ggpubr)

color_maps <- c(
  "MONTE"="#0A97B0",
  "MONTE (Mean)"="#3B7597",
  "MONTE (Weighted)"="#6FD1D7",
  "PAMES"="#F9B487",
  "PureBeta"="#A3B087",
  "InfiniumPurify"="#E8DFCA",
  "Unadjusted"="#BDA6CE"
)

cancer <- "LUSC"
metric <- "CPE"
df <- read.csv(paste0("data/benchmark/benchmark_probe_corrections_", cancer,".csv"))
df$corr_metric_group <- factor(df$corr_metric_group,
                        levels=c("[0.00, 0.25)", "[0.25, 0.50)", "[0.50, 0.75)", "[0.75, 1.00)"))
df$method <- factor(df$method,
                    levels=c("MONTE", "MONTE (Mean)", "MONTE (Weighted)", "PureBeta", "InfiniumPurify", "Unadjusted"))

# wilcox_res <- df %>%
#   group_by(corr_cpe_group) %>%
#   wilcox_test(
#     abs_cpe_pearsonr ~ method,
#     ref.group = "MONTE",
#     comparisons = list(
#       c("MONTE", "PureBeta"),
#       c("MONTE", "InfiniumPurify"),
#       c("MONTE", "Unadjusted")
#     ),
#     alternative = "less"
#   ) %>%
#   adjust_pvalue(method = "BH") %>%
#   add_significance("p.adj") %>%
#   ungroup()

# # Compute y positions for annotations
# y_pos <- df %>%
#   group_by(corr_cpe_group) %>%
#   summarise(ymax = max(abs_cpe_pearsonr, na.rm = TRUE), .groups = "drop")
# 
# wilcox_res <- wilcox_res %>%
#   left_join(y_pos, by = "corr_cpe_group") %>%
#   group_by(corr_cpe_group) %>%
#   mutate(
#     y.position = ymax + c(0.02, 0.05, 0.08)
#   ) %>%
#   ungroup()
# 
# # Optional: make labels more informative
# wilcox_res <- wilcox_res %>%
#   mutate(
#     p_label = paste0("BH-adjusted p = ", signif(p.adj, 2))
#   )

# -----------------------------
# Plot
# -----------------------------
g <- ggplot(df, aes(x = corr_metric_group, y = abs_metric_pearsonr, fill = method)) +
  geom_boxplot(outlier.size = 0, width = 0.75) +
  theme_bw() +
  labs(
    x = paste0("Unadjusted probe-purity correlation (", metric,")"),
    y = paste0("Absolute Pearson correlation\n(adjusted-probe vs. ", metric,")"),
    title = paste0("Adjusted probe-purity correlation with ", metric," (", cancer, ")"),
    fill = "Method"
  ) +
  scale_fill_manual(values = color_maps) +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )


pdf(paste0("results/supp_", cancer, "_adjusted_probe_correlation_with_purity.pdf"), width=10, height=6)
print(g)
dev.off()

