rm(list=ls())

library(ComplexHeatmap)
library(circlize)

benchmark_df <- read.csv("data/benchmark/benchmark_all_methods_purity_correlation_by_cancer.csv", check.names = FALSE)
df_test_meta <- read.csv("data/metadata/test_pan-cancer_meta.csv")
cancer_num_df <- df_test_meta %>%
  count(Cancer.type, sort = TRUE)

rownames(benchmark_df) <- benchmark_df$Cancer.type
rownames(cancer_num_df) <- cancer_num_df$Cancer.type
cancer_num_df <- cancer_num_df[order(cancer_num_df$n, decreasing = TRUE), ]
benchmark_df <- benchmark_df[cancer_num_df$Cancer.type, ]
benchmark_df$Cancer.type <- NULL
benchmark_df$n_samples <- NULL

colnames(benchmark_df)[5] <- "PAMES (pretrained)"
benchmark_df <- benchmark_df[, c("MONTE","PAMES (pretrained)","PureBeta (dataset-fitted)", "PureBeta (pretrained)","InfiniumPurify (dataset-fitted)","InfiniumPurify (pretrained)")]
benchmark_df <- benchmark_df[, c("MONTE","PAMES (pretrained)","PureBeta (pretrained)","PureBeta (dataset-fitted)","InfiniumPurify (pretrained)","InfiniumPurify (dataset-fitted)")]

mat <- t(as.matrix(benchmark_df))
col_fun <- colorRamp2(
  breaks = c(-1, -0.8, 0, 0.8, 1),
  colors = c("#18b0ca", "#c4dce4", "white", "#f99090", "#f93030")
)


n_row <- nrow(mat)
n_col <- ncol(mat)
cell_size <- unit(15, "mm")

pdf("results/fig3A_benchmark_purity_estimation.pdf", width=11, height=5)
top_annot <- HeatmapAnnotation(
  n_samples = anno_barplot(
    cancer_num_df$n,
    add_numbers = TRUE,
    numbers_gp = gpar(fontsize = 10),
    gp = gpar(fill = "#ffc3c7"),
    height = unit(2, "cm")
  ),
  annotation_label = c(n_samples = "Number of samples"),
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 12),
  annotation_name_rot = 0
)

ht <- Heatmap(
  mat,
  name = "Pearson correlation (Predicted purity vs. CPE)",
  col = col_fun,
  na_col = "gray",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  width  = unit(10, "mm") * ncol(mat),
  height = unit(10, "mm") * nrow(mat),
  
  top_annotation = top_annot,
  
  rect_gp = gpar(col = "black", lwd = 0.5),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(mat[i, j])) {
      grid.text(sprintf("%.2f", mat[i, j]), x, y,
                gp = gpar(fontsize = 10))
    }
  },
  
  heatmap_legend_param = list(
    title_position = "topcenter",
    title_gp = gpar(fontsize = 11, fontface = "plain"),
    at = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    direction = "horizontal",
    legend_width = unit(7, "cm")
  )
)

na_legend <- Legend(
  labels = "Method not supported",
  legend_gp = gpar(fill = "gray"),
  title = NULL
)

draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(na_legend),
  padding = unit(c(2, 2, 8, 2), "mm")
)

dev.off()


methods_to_compare <- setdiff(colnames(benchmark_df), "MONTE")

wilcox_results <- data.frame(
  method = character(),
  monte_mean = numeric(),
  other_mean = numeric(),
  monte_median = numeric(),
  other_median = numeric(),
  mean_difference = numeric(),
  median_difference = numeric(),
  statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (method in methods_to_compare)
{
  x <- benchmark_df$MONTE
  y <- benchmark_df[[method]]
  
  keep <- complete.cases(x, y)
  
  x <- x[keep]
  y <- y[keep]
  
  wt <- wilcox.test(
    x,
    y,
    paired = TRUE,
    alternative = "greater",
    exact = FALSE
  )
  
  wilcox_results <- rbind(
    wilcox_results,
    data.frame(
      method = method,
      monte_mean = mean(x),
      other_mean = mean(y),
      monte_median = median(x),
      other_median = median(y),
      mean_difference = mean(x - y),
      median_difference = median(x - y),
      statistic = as.numeric(wt$statistic),
      p_value = wt$p.value
    )
  )
}

# multiple testing correction
wilcox_results$FDR <- p.adjust(
  wilcox_results$p_value,
  method = "fdr"
)

# nicer formatting
wilcox_results <- wilcox_results[
  order(wilcox_results$p_value),
]