###############################################################################
#                      k-means on Lung Xenium: 3 PCs & 10 PCs
#   - Aligns barcodes across counts / PCs / spatial coords
#   - Runs k-means (K=6) for 3 PCs and 10 PCs
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

# ========================== 2) Paths & Data ==================================
work_dir   <- "/scratch/user/varogovchenko/Kmeans_runs/Lung"
data_rdata <- "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/Data_VUILD96MF/Lung_xenium_data_processed.RData"

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
out_dir <- getwd()

load(data_rdata)
# Expects:
#   counts_keep     : dgCMatrix (genes x cells), colnames = barcodes
#   pca_embeddings  : matrix/data.frame (cells x PCs), rownames = barcodes
#   spatial_keep    : tibble with columns: cell_id, x_centroid, y_centroid, Annotation_Type

# ========================== 3) Sanity & Alignment =============================
xy_mat  <- as.matrix(spatial_keep[, c("x_centroid","y_centroid")])
truth_vec <- spatial_keep$Annotation_Type

cat("Counts dim (genes x cells):", paste(dim(counts_keep), collapse=" x "), "\n")
cat("PC matrix dim (cells x PCs):", paste(dim(pca_embeddings), collapse=" x "), "\n")
cat("XY matrix dim (cells x 2):", paste(dim(xy_mat), collapse=" x "), "\n")

# ========================== 4) Helpers =======================================
plot_spatial <- function(df, color_col, title_text) {
  stopifnot(color_col %in% names(df))
  df$.grp <- factor(df[[color_col]])
  ggplot(df, aes(x = x, y = y, color = .grp)) +
    geom_point(size = 0.4, alpha = 0.9) +
    coord_equal() +
    labs(title = title_text, color = color_col, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

run_kmeans_with_pcs <- function(Y, n_pcs, K, xy, truth = NULL, dataset_tag = "Lung", out_dir = ".") {
  cat(sprintf("\n>> Running k-means (K=%d) with %d PCs ...\n", K, n_pcs))
  stopifnot(ncol(Y) >= n_pcs)
  X <- as.matrix(Y[, seq_len(n_pcs), drop = FALSE])  # cells x n_pcs

  # Robust k-means settings for stability
  set.seed(42)
  km <- kmeans(X, centers = K, nstart = 20, iter.max = 1000)
  labels <- as.integer(km$cluster)

  # ARI vs truth (if available and non-degenerate)
  ari_val <- NA_real_
  if (!is.null(truth)) {
    keep <- !is.na(truth)
    if (sum(keep) >= 2L &&
        length(unique(labels[keep])) >= 2L &&
        length(unique(truth[keep])) >= 2L) {
      truth_int <- as.integer(as.factor(truth[keep]))
      ari_val   <- mclust::adjustedRandIndex(labels[keep], truth_int)
    }
  }

  out_tbl <- tibble(
    barcode = rownames(Y),
    x       = xy[, 1],
    y       = xy[, 2],
    label   = labels,
    truth   = if (is.null(truth)) NA else truth
  )

  base <- file.path(out_dir, sprintf("KMeans_%s_%dPCs", dataset_tag, n_pcs))
  saveRDS(out_tbl, paste0(base, ".rds"))

  p <- plot_spatial(out_tbl, "label",
                    sprintf("k-means (K=%d), %s, %d PCs%s",
                            K, dataset_tag, n_pcs,
                            ifelse(is.na(ari_val), "", sprintf(", ARI = %.3f", ari_val))))
  ggsave(paste0(base, "_clusters.png"), p, width = 6, height = 5, dpi = 300)

  list(tbl = out_tbl, ari = ari_val,
       rds_path = paste0(base, ".rds"),
       png_path = paste0(base, "_clusters.png"))
}

# ========================== 5) Runs: 3 PCs & 10 PCs ===========================
K <- 6

res_10 <- run_kmeans_with_pcs(pca_embeddings, n_pcs = 10, K = K, xy = xy_mat,
                              truth = truth_vec, dataset_tag = "Lung", out_dir = out_dir)

# ========================== 6) Console Summary ================================
cat("\n========== SUMMARY (ARI excludes NA) ==========\n")
cat(sprintf("10 PCs : %s\n",  ifelse(is.na(res_10$ari), "ARI = NA", sprintf("ARI = %.6f", res_10$ari))))
cat("Saved files:\n")
cat("  - ", res_10$rds_path, "\n", sep = "")

end.time_file <- Sys.time()
file_time <- end.time_file - start.time_file
print(file_time)
cat("Finished k-means runs on Lung Xenium.\n")
###############################################################################