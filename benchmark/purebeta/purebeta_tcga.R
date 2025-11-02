# ---- Setup ----
setwd("/grain/mk98/existing_methods/purebeta/")
# renv::activate()

suppressPackageStartupMessages({
  suppressMessages({
    suppressWarnings({
      library(arrow)
      library(PureBeta)
    })
  })
})

# ---- Configuration ----
train_beta_path <- "/scratch/mk98/cancer-methyl/TCGA_Methylation_450K/processed/tumor/train"
test_beta_path  <- "/scratch/mk98/cancer-methyl/TCGA_Methylation_450K/processed/tumor/test"
meta_path       <- "/scratch/mk98/cancer-methyl/TCGA_Methylation_450K/processed/metadata/all_metadata.csv"
out_root        <- "/scratch/mk98/cancer-methyl/purebeta_results/TCGA"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_path  <- file.path(out_root, timestamp)
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

seed_num <- as.integer(Sys.getenv("PUREBETA_SEED", unset = "42"))
cores    <- max(1L, as.integer(Sys.getenv("PUREBETA_CORES", unset = "8")))

# ---- Logging ----
log_path <- file.path(out_path, "purebeta_run.log")
log_msg <- function(..., logfile = log_path) {
  msg <- paste(..., sep = "")
  timestamped <- paste0(Sys.time(), " - ", msg, "\n")
  cat(msg, "\n"); cat(timestamped, file = logfile, append = TRUE)
}

# ---- Timing helper ----
time_and_log <- function(expr, label, logfile = log_path) {
  start <- Sys.time()
  timing <- system.time(res <- tryCatch(expr, error = function(e) {
    log_msg("  Error in ", label, ": ", e$message, logfile = logfile)
    NULL
  }))
  end <- Sys.time()
  elapsed <- as.numeric(difftime(end, start, units = "secs"))
  log_msg("  ", label, " completed in ", round(elapsed, 2), "s (user=", round(timing["user"], 2),
          ", sys=", round(timing["sys"], 2), ")", logfile = logfile)
  res
}

# ---- Helpers ----
load_parquet_matrix <- function(pq_file) {
  df <- read_parquet(pq_file)
  idx_col <- "__index_level_0__"
  if (!idx_col %in% colnames(df)) stop("Missing index column in: ", pq_file)
  rownames(df) <- df[[idx_col]]
  df[[idx_col]] <- NULL
  t(as.matrix(df))  # CpGs x samples
}

# ---- Metadata ----
meta <- read.csv(meta_path, stringsAsFactors = FALSE)
stopifnot("Barcode" %in% names(meta))
rownames(meta) <- meta$Barcode
meta$CPE <- suppressWarnings(as.numeric(meta$CPE))

builtin_refs <- c("BRCA", "LUAD", "LUSC")

# ===================================================================
# 1️⃣ Build per-project reference from /train/
# ===================================================================
train_files <- list.files(train_beta_path, pattern = "_beta\\.parquet$", full.names = TRUE)
log_msg("Found ", length(train_files), " train tumor beta files.")

for (file in train_files) {
  project <- toupper(gsub("_beta\\.parquet$", "", basename(file)))
  proj_dir <- file.path(out_path, project)
  dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
  log_msg("=== [TRAIN] Building reference for ", project, " ===")
  
  tumor_beta <- tryCatch(load_parquet_matrix(file), error = function(e) {
    log_msg("  Error reading train beta: ", e$message); return(NULL)
  })
  if (is.null(tumor_beta)) next
  
  sample_ids <- intersect(colnames(tumor_beta), rownames(meta))
  if (length(sample_ids) < 20) { log_msg("  Skipping ", project, " (<20 samples)."); next }
  tumor_beta <- tumor_beta[, sample_ids, drop = FALSE]
  pur <- as.numeric(meta[sample_ids, "CPE"]); names(pur) <- sample_ids
  if (sum(is.finite(pur)) < 20) { log_msg("  Too few valid CPE purities; skipping ", project); next }
  
  ref_res <- time_and_log(
    reference_regression_generator(
      beta_values     = tumor_beta,
      tumour_purities = pur,
      set_seed        = TRUE,
      seed_num        = seed_num,
      cores           = cores
    ),
    label = paste0("[TRAIN] reference_regression_generator (", project, ")")
  )
  
  if (is.null(ref_res)) next
  saveRDS(ref_res, file.path(proj_dir, "purebeta_reference_train.rds"))
  log_msg("  Saved reference for ", project)
}

# ===================================================================
# 2️⃣ Run reference-free & reference-based modes on /test/
# ===================================================================
test_files <- list.files(test_beta_path, pattern = "_beta\\.parquet$", full.names = TRUE)
log_msg("Found ", length(test_files), " test tumor beta files.")

summary_list <- list()

for (file in test_files) {
  project <- toupper(gsub("_beta\\.parquet$", "", basename(file)))
  proj_dir <- file.path(out_path, project)
  dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
  log_msg("=== [TEST] Running both modes for ", project, " ===")
  
  proj_start <- Sys.time()
  
  tumor_beta <- tryCatch(load_parquet_matrix(file), error = function(e) {
    log_msg("  Error reading test beta: ", e$message); return(NULL)
  })
  if (is.null(tumor_beta)) next
  
  sample_ids <- intersect(colnames(tumor_beta), rownames(meta))
  if (length(sample_ids) < 5) { log_msg("  Skipping ", project, " (<5 samples)."); next }
  tumor_beta <- tumor_beta[, sample_ids, drop = FALSE]
  meta_sub <- meta[sample_ids, , drop = FALSE]
  true_pur <- meta_sub$CPE; names(true_pur) <- rownames(meta_sub)
  
  # ---- Reference-free ----
  ref_free <- time_and_log(
    beta_correction_for_cohorts(
      beta_values = tumor_beta,
      set_seed    = TRUE,
      seed_num    = seed_num,
      cores       = cores
    ),
    label = paste0("[TEST] ref-free beta_correction_for_cohorts (", project, ")")
  )
  
  purity_free <- if (!is.null(ref_free)) ref_free$estimated_purities else NULL
  adj_free    <- if (!is.null(ref_free)) ref_free$adjusted_betas else NULL
  if (!is.null(ref_free)) saveRDS(ref_free, file.path(proj_dir, "purebeta_refFree_results.rds"))
  
  # ---- Reference-based ----
  ref_based <- NULL; purity_based <- NULL; adj_based <- NULL
  use_builtin <- project %in% builtin_refs
  ref_path <- file.path(proj_dir, "purebeta_reference_train.rds")
  has_train_ref <- file.exists(ref_path)
  
  if (use_builtin || has_train_ref) {
    ref_model <- if (has_train_ref) readRDS(ref_path) else NULL
    pur_res <- time_and_log(
      if (use_builtin && is.null(ref_model)) {
        purity_estimation(beta_values = tumor_beta, tumour_type = project, cores = cores)
      } else {
        purity_estimation(reference_regressions = ref_model, beta_values = tumor_beta, cores = cores)
      },
      label = paste0("[TEST] reference-based purity_estimation (", project, ")")
    )
    
    if (!is.null(pur_res)) {
      purity_based <- pur_res$estimated_purities
      saveRDS(pur_res, file.path(proj_dir, "purebeta_refBased_purity.rds"))
    }
    
    if (!is.null(purity_based)) {
      adj_based <- time_and_log(
        reference_based_beta_correction(
          beta_values           = tumor_beta,
          purities              = pur_res,
          reference_regressions = if (!is.null(ref_model)) ref_model else NULL,
          refitting             = FALSE,
          cores                 = cores
        ),
        label = paste0("[TEST] reference_based_beta_correction (", project, ")")
      )
      if (!is.null(adj_based))
        saveRDS(adj_based, file.path(proj_dir, "purebeta_refBased_adjustedBetas.rds"))
    }
  } else {
    log_msg("  No built-in or trained reference for ", project, "; skipping ref-based.")
  }
  
  # ---- Evaluation ----
  calc_metrics <- function(est, true) {
    if (is.null(est) || all(is.na(est))) return(c(cor = NA_real_, mse = NA_real_))
    common <- intersect(names(est), names(true))
    if (length(common) < 3) return(c(cor = NA_real_, mse = NA_real_))
    est <- est[common]; tru <- true[common]
    c(cor = suppressWarnings(cor(est, tru, use="complete.obs")),
      mse = mean((est - tru)^2, na.rm=TRUE))
  }
  
  m_free  <- calc_metrics(purity_free,  true_pur)
  m_based <- calc_metrics(purity_based, true_pur)
  
  proj_end <- Sys.time()
  proj_elapsed <- as.numeric(difftime(proj_end, proj_start, units = "secs"))
  
  summary_list[[project]] <- data.frame(
    project = project,
    n_samples = ncol(tumor_beta),
    mean_refFree  = mean(purity_free,  na.rm=TRUE),
    cor_refFree   = m_free[["cor"]],
    mse_refFree   = m_free[["mse"]],
    mean_refBased = mean(purity_based, na.rm=TRUE),
    cor_refBased  = m_based[["cor"]],
    mse_refBased  = m_based[["mse"]],
    runtime_sec   = round(proj_elapsed, 2),
    stringsAsFactors = FALSE
  )
  
  log_msg("  Completed ", project, " (", round(proj_elapsed, 2), "s): mean_refFree=",
          round(mean(purity_free,na.rm=TRUE),3),
          ", mean_refBased=", round(mean(purity_based,na.rm=TRUE),3))
}

# ===================================================================
# 3️⃣ Save summary
# ===================================================================
summary_df <- if (length(summary_list)) do.call(rbind, summary_list) else data.frame()
write.csv(summary_df, file.path(out_path, "purebeta_summary.csv"), row.names = FALSE)
log_msg("Saved summary to ", out_path)
log_msg("All projects completed.")
