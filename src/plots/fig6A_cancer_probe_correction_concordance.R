rm(list=ls())
library(tidyverse)
library(RColorBrewer)

library(ggpubr)
library(ggrepel)

df <- read_csv("data/application/SUMMARY_all_cancers_splithalf_tumor_only.csv") %>%
  filter(!cancer %in% c("GBM", "CESC", "SKCM")) %>%
  arrange(split_pearson_median_uncorrected) %>%
  mutate(cancer = factor(cancer, levels = cancer))

cancer_colors <- colorRampPalette(
  brewer.pal(12, "Paired")
)(n_distinct(df$cancer))

names(cancer_colors) <- levels(df$cancer)
cancer_colors["BRCA"] <- "#FFC94D"

df_long <- df %>%
  select(
    cancer,
    split_pearson_median_uncorrected,
    split_pearson_median_corrected
  ) %>%
  pivot_longer(
    cols = c(
      split_pearson_median_uncorrected,
      split_pearson_median_corrected
    ),
    names_to = "type",
    values_to = "split_r_median"
  ) %>%
  mutate(
    type = recode(
      type,
      split_pearson_median_uncorrected = "Uncorrected",
      split_pearson_median_corrected = "Corrected"
    ),
    type = factor(
      type,
      levels = c("Uncorrected", "Corrected")
    )
  )

p <- ggplot(
  df_long,
  aes(
    x = type,
    y = split_r_median,
    group = cancer
  )
) +
  geom_line(
    aes(color = cancer),
    linewidth = 0.8
  ) +
  geom_point(
    aes(color = cancer),
    size = 2.5
  ) +
  geom_text_repel(
    data = df_long %>%
      filter(type == "Corrected"),
    aes(
      label = cancer,
      color = cancer
    ),
    nudge_x = 0.3,
    direction = "y",
    hjust = 0,
    size = 5,
    segment.color = "grey60",
    segment.size = 0.3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  stat_compare_means(
    comparisons = list(
      c("Uncorrected", "Corrected")
    ),
    paired = TRUE,
    method = "wilcox.test",
    alternative = "less",
    label = "p.signif",
    size=
  ) +
  scale_color_manual(
    values = cancer_colors
  ) +
  labs(
    x = NULL,
    y = expression(
      "Concordance (Pearson " * rho * ")"
    ),
    color = "Cancer"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    plot.margin = margin(
      5.5,
      80,
      5.5,
      5.5
    ),
    plot.title = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "NONE",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.background = element_blank()
  ) +
  guides(
    color = guide_legend(nrow = 4)
  )

pdf("results/fig6A_cancer_probe_correction_concordance.pdf", width=6, height=7)
p
dev.off()
