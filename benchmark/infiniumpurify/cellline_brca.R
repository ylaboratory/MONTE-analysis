setwd("/grain/mk98/existing_methods/infiniumpurify")

# ---------- Setup ----------
library(arrow)      # for read_parquet()
library(readr)      # for read_csv()
library(dplyr)
library(tibble)
library(ggplot2)

cancer <- "BRCA"

train_path <- "/grain/mk98/cancer-methyl/TCGA_Methylation_450K/processed/tumor/train"
test_path  <- "/grain/mk98/cancer-methyl/TCGA_Methylation_450K/processed/tumor/test"
meta_path  <- "/grain/mk98/cancer-methyl/TCGA_Methylation_450K/processed/metadata/tumor_metadata.csv"
xy_file    <- "/grain/mk98/cancer-methyl/probe_selection_files/xy_probes.txt"

# ---------- Load metadata and XY probes ----------
xy_probes <- read_csv(xy_file, col_names = FALSE, show_col_types = FALSE)[[1]]
all_metadata <- read_csv(meta_path, show_col_types = FALSE) %>%
  column_to_rownames(var = colnames(.)[1])

# ---------- Load beta matrices ----------
idx_col <- "__index_level_0__"

train_beta <- read_parquet(file.path(train_path, paste0(cancer, "_beta.parquet"))) %>%
  as.data.frame()                         # convert tibble -> data.frame first
rownames(train_beta) <- train_beta[[idx_col]]
train_beta[[idx_col]] <- NULL

test_beta <- read_parquet(file.path(test_path, paste0(cancer, "_beta.parquet"))) %>%
  as.data.frame()
rownames(test_beta) <- test_beta[[idx_col]]
test_beta[[idx_col]] <- NULL

# Now both train_beta and test_beta are plain data.frames with numeric columns
cat("train_beta dim:", dim(train_beta), "\n")
cat("test_beta dim:",  dim(test_beta),  "\n")
cat("example rownames:", head(rownames(test_beta)), "\n")

# Align metadata
train_meta <- all_metadata[rownames(train_beta), , drop = FALSE]
test_meta  <- all_metadata[rownames(test_beta),  , drop = FALSE]

# ---------- Compute M-values ----------
beta_to_m <- function(beta) log2(beta / (1 - beta))
train_m <- beta_to_m(train_beta)
test_m  <- beta_to_m(test_beta)

cat("Original shapes:\n")
cat("train_m:", dim(train_m), "\n")
cat("test_m:",  dim(test_m),  "\n")

# ---------- Remove NA samples and XY probes ----------
train_m <- train_m[complete.cases(train_m), !(colnames(train_m) %in% xy_probes)]
test_m  <- test_m[complete.cases(test_m),  !(colnames(test_m)  %in% xy_probes)]

cat("Filtered shapes:\n")
cat("train_m:", dim(train_m), "\n")
cat("test_m:",  dim(test_m),  "\n")

# ---------- Load cell line data ----------
cellline_path <- "/grain/mk98/cancer-methyl/external/celline_GSE68379"

cellline_beta <- read_parquet(file.path(cellline_path, "GSE68379_all_beta.parquet")) %>%
  as.data.frame()
rownames(cellline_beta) <- cellline_beta$probe_id
cellline_beta$probe_id <- NULL
cellline_beta <- t(cellline_beta)

cellline_meta <- read_parquet(file.path(cellline_path, "GSE68379_all_meta.parquet")) %>%
  as.data.frame()
rownames(cellline_meta) <- cellline_meta$sample_id

# Align meta rows to beta rows
cellline_meta <- cellline_meta[rownames(cellline_beta), , drop = FALSE]

# Remove XY probes and NA columns
cellline_beta <- cellline_beta[, !(colnames(cellline_beta) %in% xy_probes)]
cellline_beta <- cellline_beta[, colSums(is.na(cellline_beta)) == 0]

# Compute M-values
cellline_m <- beta_to_m(cellline_beta)

# Filter for chosen cancer type (BRCA)
cellline_meta_cancer <- subset(cellline_meta, project == cancer)
cellline_m_cancer <- cellline_m[rownames(cellline_meta_cancer), , drop = FALSE]

cat("Cell line BRCA shapes:\n")
cat("cellline_m_cancer:", dim(cellline_m_cancer), "\n")
cat("cellline_meta_cancer:", dim(cellline_meta_cancer), "\n")

library(InfiniumPurify)

cellline_vec <- getPurity(tumor.data = as.data.frame(t(cellline_beta_cancer)),
                          normal.data = NULL,
                          tumor.type = "BRCA")

ggplot(data.frame(purity = cellline_vec), aes(x = purity)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  # geom_vline(aes(xintercept = mean(purity, na.rm = TRUE)),
  #            color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Tumor Purity Estimates (BRCA): cell line",
       x = "Estimated Purity",
       y = "Number of Samples") +
  theme_minimal(base_size = 13)

tcga_vec <- getPurity(tumor.data = as.data.frame(t(test_beta)),
                      normal.data = NULL,
                      tumor.type = "BRCA")
df <- data.frame(True = test_meta$CPE, Pred = tcga_vec)
df <- df[complete.cases(df), ]
corr <- cor(df$True, df$Pred)
mse  <- mean((df$True - df$Pred)^2)

# Scatter plot
ggplot(df, aes(x = True, y = Pred)) +
  geom_point(alpha = 0.5, color = "royalblue") +
  geom_abline(linetype = "dashed") +
  labs(
    title = sprintf("%s: corr=%.3f, mse=%.4f", cancer, corr, mse),
    x = "True CPE",
    y = "Predicted Purity"
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal()

