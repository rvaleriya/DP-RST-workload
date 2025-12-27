###############################################################################
#                  SC.MEB on Lung Xenium: 3 PCs & 10 PCs (K = 6)
#   - Aligns barcodes across PCs and spatial coords
#   - Builds spatial neighborhood (radius cutoff = 150)
#   - Runs SC.MEB for 3 PCs and 10 PCs with K=6
#   - Computes ARI vs Annotation_Type (if present), saves RDS + PNGs
###############################################################################

writeLines("R script started", stderr())
start.time_file <- Sys.time()
set.seed(42)

# ========================== 1) Libraries & LibPath ============================
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(Matrix)
library(dplyr)
library(ggplot2)
library(mclust)
library(tibble)
library(SC.MEB)
library(FNN)
library(parallel)

# ========================== 2) Paths & Data ==================================
work_dir   <- "/scratch/user/varogovchenko/SCMEB_runs/Lung"
data_rdata <- "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/Data_VUILD96MF/Lung_xenium_data_processed.RData"

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
out_dir <- getwd()

load(data_rdata)
# Expects:
#   counts_keep     : dgCMatrix (genes x cells), colnames = barcodes
#   pca_embeddings  : matrix/data.frame (cells x PCs), rownames = barcodes
#   spatial_keep    : tibble/data.frame with columns:
#                     cell_id, x_centroid, y_centroid, Annotation_Type

stopifnot(exists("pca_embeddings"), exists("spatial_keep"))

# ========================== 3) Sanity & Alignment =============================
# Extract XY and truth
xy_mat    <- as.matrix(spatial_keep[, c("x_centroid","y_centroid")])
truth_vec <- spatial_keep$Annotation_Type
cat("PC matrix dim (cells x PCs):", paste(dim(pca_embeddings), collapse=" x "), "\n")
cat("XY matrix dim (cells x 2):",   paste(dim(xy_mat), collapse=" x "), "\n")

# ========================== 4) SC.MEB Hyper-parameters ========================
beta_grid     <- seq(0, 5, 0.2)
K_set         <- 6       
parallel_run  <- TRUE
num_core      <- max(1L, min(8L, detectCores()))
PX            <- TRUE
maxIter_ICM   <- 10
maxIter       <- 50

# ========================== 5) Spatial Neighborhood ===========================
# Xenium spatial scale: start with 150; adjust if too sparse/dense.
Adj_sp <- getneighborhood_fast(xy_mat, cutoff = 150)

# ========================== 6) Helpers =======================================
extract_labels <- function(scmeb_fit_singleK) {
  # SC.MEB returns an object whose $x holds labels; handle vector/matrix cases
  x_obj <- scmeb_fit_singleK$x
  if (is.matrix(x_obj)) as.integer(x_obj[, 1]) else as.integer(x_obj)
}

plot_spatial <- function(df, color_col, title_text) {
  stopifnot(color_col %in% names(df))
  df$.grp <- factor(df[[color_col]])
  ggplot(df, aes(x = x, y = y, color = .grp)) +
    geom_point(size = 0.4, alpha = 0.9) +
    coord_fixed() +
    labs(title = title_text, color = color_col, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

run_scmeb_once <- function(pcs_mat, n_pcs, K_val, xy, truth = NULL,
                           dataset_tag = "Lung", out_dir = ".") {
  message(sprintf(">> Running SC.MEB with %d PCs, K = %d ...", n_pcs, K_val))

  stopifnot(ncol(pcs_mat) >= n_pcs)
  X <- as.matrix(pcs_mat[, seq_len(n_pcs), drop = FALSE])  # cells x n_pcs

  fit <- SC.MEB(
    y            = X,
    Adj_sp       = Adj_sp,
    beta_grid    = beta_grid,
    K_set        = K_set,          # length 1: K_val
    parallel     = parallel_run,
    num_core     = num_core,
    PX           = PX,
    maxIter_ICM  = maxIter_ICM,
    maxIter      = maxIter
  )

  # Select result for K=K_val (first column since K_set has length 1)
  resK   <- fit[, 1]
  labels <- extract_labels(resK)

  # Compute ARI vs truth (if feasible)
  ari_val <- NA_real_
  if (!is.null(truth)) {
    keep <- !is.na(truth)
    if (sum(keep) >= 2L &&
        length(unique(labels[keep])) >= 2L &&
        length(unique(truth[keep])) >= 2L) {
      ref <- as.integer(as.factor(truth[keep]))
      ari_val <- adjustedRandIndex(labels[keep], ref)
    }
  }

  out_tbl <- tibble(
    barcode = rownames(pcs_mat),
    x       = xy[, 1],
    y       = xy[, 2],
    label   = labels,
    truth   = if (is.null(truth)) NA else truth
  )

  base <- file.path(out_dir, sprintf("SC-MEB_%s_%dPCs_K%d", dataset_tag, n_pcs, K_val))
  rds_path <- paste0(base, ".rds")
  png_path <- paste0(base, "_clusters.png")

  saveRDS(out_tbl, rds_path)

  title_txt <- sprintf("SC.MEB (K=%d), %s, %d PCs%s",
                       K_val, dataset_tag, n_pcs,
                       ifelse(is.na(ari_val), "", sprintf(", ARI = %.3f", ari_val)))
  p <- plot_spatial(out_tbl, "label", title_txt)
  ggsave(png_path, p, width = 6, height = 5, dpi = 300)

  list(tbl = out_tbl, ari = ari_val, rds_path = rds_path, png_path = png_path)
}

# ========================== 7) Runs: 3 PCs & 10 PCs ===========================
K <- 6
res_10 <- run_scmeb_once(pca_embeddings, n_pcs = 10, K_val = K, xy = xy_mat,
                         truth = truth_vec, dataset_tag = "Lung", out_dir = out_dir)

# ========================== 8) Console Summary ================================
cat("\n========== SUMMARY (ARI excludes NA) ==========\n")
cat(sprintf("10 PCs : %s\n",  ifelse(is.na(res_10$ari), "ARI = NA", sprintf("ARI = %.6f", res_10$ari))))
cat("Saved files:\n")
cat("  - ", res_10$rds_path, "\n", sep = "")

end.time_file <- Sys.time()
file_time <- end.time_file - start.time_file
print(file_time)
cat("Finished SC.MEB runs on Lung Xenium.\n")
###############################################################################