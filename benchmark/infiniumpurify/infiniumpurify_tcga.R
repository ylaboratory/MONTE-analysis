setwd('/grain/mk98/existing_methods/infiniumpurify/')
# renv::activate()

suppressPackageStartupMessages({
  suppressMessages({
    suppressWarnings({
      library(arrow)
      library(InfiniumPurify)
    })
  })
})

plot_probe_correlation <- FALSE  # set TRUE to compute/plot probe–purity correlations
train_tumor_path <- "/scratch/mk98/cancer-methyl/TCGA_Methylation_450K/processed/tumor/train"
test_tumor_path <- "/scratch/mk98/cancer-methyl/TCGA_Methylation_450K/processed/tumor/test"
train_normal_path <- "/scratch/mk98/cancer-methyl/TCGA_Methylation_450K/processed/normal/train"
test_normal_path <- "/scratch/mk98/cancer-methyl/TCGA_Methylation_450K/processed/normal/test"

all_meta_path <- "/scratch/mk98/cancer-methyl/TCGA_Methylation_450K/processed/metadata/all_metadata.csv"
out_path_root <- "/scratch/mk98/cancer-methyl/infiniumpurify_results/TCGA"
xy_probes_df <- read.csv("/scratch/mk98/cancer-methyl/probe_selection_files/xy_probes.txt", header = FALSE, stringsAsFactors = FALSE)
xy_probes <- xy_probes_df[[1]]

load_beta_parquet <- function(path) {
  df <- read_parquet(path)
  rownames(df) <- df$`__index_level_0__`
  df$`__index_level_0__` <- NULL
  mat <- t(as.matrix(df))
  mat
}

filter_align_cpgs <- function(
  tumor_beta, 
  normal_beta = NULL, 
  xy_probes = character(),
  project = "<project>"
) {
  tumor_keep <- (complete.cases(tumor_beta) | !rownames(tumor_beta) %in% xy_probes)
  if (!is.null(normal_beta)) {
    normal_keep <- (complete.cases(normal_beta) | !rownames(normal_beta) %in% xy_probes)
    valid_cpgs <- intersect(rownames(tumor_beta)[tumor_keep], rownames(normal_beta)[normal_keep])
    tumor_beta  <- tumor_beta[valid_cpgs, , drop = FALSE]
    normal_beta <- normal_beta[valid_cpgs, , drop = FALSE]
  } else {
    valid_cpgs <- rownames(tumor_beta)[tumor_keep]
    tumor_beta <- tumor_beta[valid_cpgs, , drop = FALSE]
    normal_beta <- NULL
  }
  list(skipped = FALSE, valid_cpgs = valid_cpgs,
       tumor_beta = tumor_beta, normal_beta = normal_beta)
}

calc_metrics <- function(est, tru) {
  if (is.null(est)) return(c(cor = NA_real_, mse = NA_real_))
  common <- intersect(names(est), names(tru))
  if (length(common) < 2) return(c(cor = NA_real_, mse = NA_real_))
  est <- est[common]
  tru <- tru[common]
  # Require at least two complete (non-NA, finite) pairs for Pearson correlation
  c(
    cor = cor(est, tru, method = "pearson"),
    mse = mean((est - tru)^2, na.rm = TRUE)
  )
}
# optional mapping for tumor types not in iDMC (e.g., READ)
tumor_type_map <- c(READ = "COAD")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path <- file.path(out_path_root, timestamp)
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

# --- Logging ---
log_path <- file.path(out_path, "infiniumpurify_run.log")
log_msg <- function(..., logfile = log_path) {
  msg <- paste(..., sep = "")
  timestamped <- paste0(Sys.time(), " - ", msg, "\n")
  cat(msg, "\n")
  cat(timestamped, file = logfile, append = TRUE)
}

# --- Load metadata ---
meta <- read.csv(all_meta_path, stringsAsFactors = FALSE)
stopifnot("Barcode" %in% names(meta))
rownames(meta) <- meta$Barcode
meta$CPE <- suppressWarnings(as.numeric(meta$CPE))

# --- Discover tumor files ---
train_tumor_files <- list.files(train_tumor_path, pattern = "_beta\\.parquet$", full.names = TRUE)
log_msg("Found ", length(train_tumor_files), " train tumor files.", logfile = log_path)
test_tumor_files <- list.files(test_tumor_path, pattern = "_beta\\.parquet$", full.names = TRUE)
log_msg("Found ", length(test_tumor_files), " test tumor files.", logfile = log_path)

train_normal_files <- list.files(train_normal_path, pattern = "_beta\\.parquet$", full.names = TRUE)
log_msg("Found ", length(train_normal_files), " train normal files.", logfile = log_path)
test_normal_files <- list.files(test_normal_path, pattern = "_beta\\.parquet$", full.names = TRUE)
log_msg("Found ", length(test_normal_files), " test normal files.", logfile = log_path)

summary_list <- list()

# Preload iDMC names for safe tumor.type mapping (if available)
supp_types <- tryCatch({
  data(iDMC, package = "InfiniumPurify")
  names(iDMC)
}, error = function(e) character(0))

for (file in test_tumor_files) {
  project <- gsub("_beta\\.parquet$", "", basename(file))
  log_msg("=== Processing ", project, " ===", logfile = log_path)
  
  # --- Load tumor beta (samples x probes in parquet; transpose to CpGs x samples) ---
  test_tumor_beta <- load_beta_parquet(file)

  # --- Load normal beta if exists ---
  normal_file <- file.path(test_normal_path, paste0(project, "_beta.parquet"))
  test_normal_beta <- NULL
  if (file.exists(normal_file)) {
    test_normal_beta <- load_beta_parquet(normal_file)
  }

  # load train tumor and normal beta and meta
  train_tumor_file <- file.path(train_tumor_path, paste0(project, "_beta.parquet"))
  train_tumor_beta <- NULL
  if (file.exists(train_tumor_file)) {
    train_tumor_beta <- load_beta_parquet(train_tumor_file)
  }
  train_normal_file <- file.path(train_normal_path, paste0(project, "_beta.parquet"))
  train_normal_beta <- NULL
  if (file.exists(train_normal_file)) {
    train_normal_beta <- load_beta_parquet(train_normal_file)
  }

  # Align TRAIN tumor/normal together to define the reference CpG set
  res <- filter_align_cpgs(
    train_tumor_beta, 
    if (exists("train_normal_beta")) train_normal_beta else NULL,
    xy_probes = xy_probes,
    project = project
  )

  if (res$skipped) next
  train_tumor_beta  <- res$tumor_beta
  train_normal_beta <- res$normal_beta
  valid_cpgs_train <- res$valid_cpgs

  res <- NULL

  # Align TEST matrices to the TRAIN-defined CpG set (no recomputation of CpGs from test)
  valid_cpgs_test <- intersect(valid_cpgs_train, rownames(test_tumor_beta))
  test_tumor_beta <- test_tumor_beta[valid_cpgs_test, , drop = FALSE]
  if (exists("test_normal_beta") && !is.null(test_normal_beta)) {
    test_normal_beta <- test_normal_beta[intersect(valid_cpgs_test, rownames(test_normal_beta)), , drop = FALSE]
  }

  train_tumor_meta <- meta[intersect(rownames(meta), colnames(train_tumor_beta)), , drop = FALSE]
  test_tumor_meta  <- meta[intersect(rownames(meta), colnames(test_tumor_beta)), , drop = FALSE]
  
  cat(dim(train_tumor_beta), dim(train_normal_beta))
  cat(dim(test_tumor_beta), dim(test_normal_beta))

  # --- Output dir ---
  proj_out_dir <- file.path(out_path, project)
  if (!dir.exists(proj_out_dir)) dir.create(proj_out_dir, recursive = TRUE, showWarnings = FALSE)

  # train our own reference + estimate purity
  purity_rf <- NULL
  purity_rb <- NULL
  if (!is.null(train_normal_beta) && ncol(train_normal_beta) >= 20 && ncol(train_tumor_beta) >= 20) {
    log_msg("Running reference-free getPurity() for ", project, logfile = log_path)
    common_cpgs <- intersect(rownames(train_normal_beta), rownames(test_tumor_beta))
    train_normal_beta_rf <- train_normal_beta[common_cpgs, , drop = FALSE]
    train_tumor_beta_rf  <- train_tumor_beta[common_cpgs, , drop = FALSE]

    purity_rf <- InfiniumPurify::getPurity(tumor.data = train_tumor_beta_rf, normal.data = train_normal_beta_rf)
    if (!is.null(purity_rf)) saveRDS(purity_rf, file.path(proj_out_dir, "infiniumpurify_refFree_purity.rds"))
  } else {
    log_msg("Skipping ref-free for ", project, ": need >=20 tumors and normals in train.", logfile = log_path)
  }
  
  # use reference + estimate purity
  tumor_label <- toupper(project)
  if (length(supp_types)) {
    if (!tumor_label %in% supp_types && tumor_label %in% names(tumor_type_map)) {
      log_msg("Mapping tumor.type ", tumor_label, " → ", tumor_type_map[[tumor_label]], " for iDMC.", logfile = log_path)
      tumor_label <- tumor_type_map[[tumor_label]]
    }
  }
  log_msg("Running reference-based getPurity() for ", project, " (type=", tumor_label, ")", logfile = log_path)
  purity_rb <- InfiniumPurify::getPurity(tumor.data = test_tumor_beta, tumor.type = tumor_label)
  if (!is.null(purity_rb)) saveRDS(purity_rb, file.path(proj_out_dir, "infiniumpurify_refBased_purity.rds"))

  mean_free  <- if (!is.null(purity_rf)) mean(purity_rf,  na.rm = TRUE) else NA_real_
  mean_based <- if (!is.null(purity_rb)) mean(purity_rb, na.rm = TRUE) else NA_real_
  log_msg("Mean purity ref-free=", round(mean_free, 3), " ; ref-based=", round(mean_based, 3), logfile = log_path)
  
  # --- Metrics vs CPE (optional) ---
  true_pur <- as.numeric(test_tumor_meta$CPE)
  names(true_pur) <- rownames(test_tumor_meta)
  
  train_true_pur <- as.numeric(train_tumor_meta$CPE)
  names(train_true_pur) <- rownames(train_tumor_meta)
  
  res_free  <- if (!is.null(purity_rf)) calc_metrics(purity_rf, train_true_pur) else c(cor = NA_real_, mse = NA_real_)
  res_based <- if (!is.null(purity_rb)) calc_metrics(purity_rb, true_pur) else c(cor = NA_real_, mse = NA_real_)
  log_msg("Ref-free purity: cor=", round(res_free[["cor"]], 3), " ; mse=", round(res_free[["mse"]], 3), logfile = log_path)
  log_msg("Ref-based purity: cor=", round(res_based[["cor"]], 3), " ; mse=", round(res_based[["mse"]], 3), logfile = log_path)

  # --- Beta adjustment (REQUIRES normal.data) ---
  adj_free <- adj_based <- NULL
  if (!is.null(purity_rf) && !is.null(test_normal_beta)) {
    adj_free <-tryCatch({
      InfiniumPurify::InfiniumPurify(tumor.data = test_tumor_beta, normal.data = test_normal_beta, purity = purity_rf)
    }, error = function(e) { log_msg("Adj ref-free error: ", e$message, logfile = log_path); NULL })
    if (!is.null(adj_free))
      saveRDS(adj_free, file.path(proj_out_dir, "infiniumpurify_refFree_adjBeta.rds"))
  } else if (!is.null(purity_rf) && is.null(test_normal_beta)) {
    log_msg("Ref-free purity available but no normals for adjustment; skipping adj (", project, ").", logfile = log_path)
  }
  
  if (!is.null(purity_rb)) {
    if (!is.null(test_normal_beta)) {
      adj_based <- tryCatch({
        InfiniumPurify::InfiniumPurify(tumor.data = test_tumor_beta, normal.data = test_normal_beta, purity = purity_rb)
      }, error = function(e) { log_msg("Adj ref-based error: ", e$message, logfile = log_path); NULL })
      if (!is.null(adj_based))
        saveRDS(adj_based, file.path(proj_out_dir, "infiniumpurify_refBased_adjBeta.rds"))
    } else {
      log_msg("No normal data available; skipping reference-based adjustment for ", project, ".", logfile = log_path)
    }
  }
  
  # --- Summary row ---
  summary_list[[project]] <- data.frame(
    project        = project,
    n_tumor_train  = ncol(train_tumor_beta),
    n_normal_train = ifelse(is.null(train_normal_beta), 0L, ncol(train_normal_beta)),
    n_tumor_test   = ncol(test_tumor_beta),
    n_normal_test  = ifelse(is.null(test_normal_beta), 0L, ncol(test_normal_beta)),
    mean_ref_free  = mean_free,
    cor_ref_free   = res_free[["cor"]],
    mse_ref_free   = res_free[["mse"]],
    mean_ref_based = mean_based,
    cor_ref_based  = res_based[["cor"]],
    mse_ref_based  = res_based[["mse"]],
    stringsAsFactors = FALSE
  )
  
  log_msg("Completed ", project, ".", logfile = log_path)
}

# --- Save summary ---
summary_df <- if (length(summary_list)) do.call(rbind, summary_list) else data.frame()
write.csv(summary_df, file.path(out_path, "infiniumpurify_summary.csv"), row.names = FALSE)
log_msg("Saved summary to ", out_path, logfile = log_path)
log_msg("All projects completed.", logfile = log_path)

cat("reference-free uses training set to measure performance due to model structure")
summary_df
