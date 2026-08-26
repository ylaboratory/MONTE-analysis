rm(list=ls())

library(ggplot2)
library(ggthemes)
library(tidyverse)

cancer <- "BRCA"
df <- read.csv(paste0("data/monte_outputs/pancancer/monte_", cancer,"_corrected_probe_correlations.csv"))
colnames(df) <- c("CpG_ID", "Uncorrected", "MONTE corrected")
color_map <- c("MONTE corrected"="#0A97B0", "Uncorrected"="#BDA6CE")
df_long <- df %>%
  select(Uncorrected, `MONTE corrected`) %>%
  pivot_longer(cols = everything(),
               names_to = "Type",
               values_to = "Value")

g_density <- ggplot(df_long, aes(x = Value, color = Type, fill = Type)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  theme_bw() +
  labs(x = "Pearson correlation", y = "Density", 
       title=paste0("Probe correlation with purity (", cancer, ")"), fill="", color="") +
  scale_color_manual(values=color_map) +
  scale_fill_manual(values=color_map) +
  theme(
    plot.title = element_text(size = 16, hjust=0.5),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    legend.position = c(.15, .85),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

g_hist <- ggplot(df_long, aes(x = Value, fill = Type)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 50) +
  theme_bw() +
  labs(x = "Pearson correlation (Probe value vs. CPE)", y = "Frequency",
       title=paste0("Probe correlation with CPE (", cancer, ")"), fill="", color="") +
  scale_color_manual(values=color_map) +
  scale_fill_manual(values=color_map) +
  theme(
    plot.title = element_text(size = 18, hjust=0.5),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    # legend.background = element_rect(fill = "transparent", color = NA),
    # legend.box.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

pdf(paste0("results/fig5A_adjusted_probes_correlation_purity_", cancer,".pdf"), width=5.5, height=5)
print(g_hist)
dev.off()
