# ========================== BASS on Lung Xenium (10 PCs) ==========================
# Reproducible, end-to-end run: align barcodes, run BASS, save outputs, print ARI.
# ================================================================================

writeLines("R script started", stderr())
start.time_file <- Sys.time()
set.seed(42)

# Set a writable library path
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

# Libraries
library(ggplot2)
library(dplyr)
library(mclust)
library(BASS)
library(Seurat)
library(readr)
library(Matrix)

# Paths
setwd("/scratch/user/varogovchenko/BASS_runs/Lung")
load("/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/Data_VUILD96MF/Lung_xenium_data_processed.RData")
# Expects objects: counts_keep (dgCMatrix genes x cells), pca_embeddings (cells x PCs), spatial_keep (tibble)

# ---- Align barcodes across counts / PCs / spatial coords -----------------------
barcodes_counts <- colnames(counts_keep)
barcodes_pcs    <- rownames(pca_embeddings)
barcodes_spat   <- spatial_keep$cell_id
barcodes <- Reduce(intersect, list(barcodes_counts, barcodes_pcs, barcodes_spat))

print("Counts matrix dimension")
print(dim(counts_keep))

print("Barcodes dimension")
print(dim(barcodes))

count_t <- counts_keep[, barcodes, drop = FALSE]                            # genes x cells
Y       <- pca_embeddings[barcodes, , drop = FALSE]                         # cells x PCs
xy_mat  <- as.matrix(spatial_keep[match(barcodes, spatial_keep$cell_id),
                                  c("x_centroid","y_centroid")])
rownames(xy_mat) <- barcodes
z <- spatial_keep$Annotation_Type[match(barcodes, spatial_keep$cell_id)]

# ---- Helper: quick scatter for clusters/labels --------------------------------
plot_spatial <- function(df, title_text) {
  ggplot(df, aes(x = x, y = y, color = factor(z))) +
    geom_point(size = 0.4) +
    coord_equal() +
    labs(title = title_text, color = "Domain") +
    theme_minimal()
}

# ---- Run BASS with N PCs -------------------------------------------------------
run_bass_with_pcs <- function(Y, n_pcs, count_t, xy_mat, z_vec = NULL, dataset_tag = "Lung") {
  bass <- createBASSObject(
    X  = list(count_t),        # genes x cells
    xy = list(xy_mat),         # cells x 2 (rownames = barcodes)
    C  = 20, 
    R  = 6,
    beta_method = "SW"
  )

  pc_mat <- t(as.matrix(Y[, seq_len(n_pcs), drop = FALSE]))  # PCs x cells
  colnames(pc_mat) <- colnames(count_t)
  bass@X_run <- pc_mat

  bass <- BASS.run(bass)
  bass <- BASS.postprocess(bass, adjustLS = TRUE)

  zlabels <- as.integer(bass@results$z[[1]])

  df_bass <- tibble(
    barcode = colnames(count_t),
    x = xy_mat[, 1],
    y = xy_mat[, 2],
    z = zlabels
  )

  ari_val <- NA_real_
  if (!is.null(z_vec)) {
    keep <- !is.na(z_vec)
    if (sum(keep) >= 2L && length(unique(zlabels[keep])) >= 2L && length(unique(z_vec[keep])) >= 2L) {
      truth_int <- as.integer(as.factor(z_vec[keep]))
      ari_val   <- mclust::adjustedRandIndex(zlabels[keep], truth_int)
    }
  }

  base <- paste0("BASS_", dataset_tag, "_", n_pcs, "PCs")
  saveRDS(df_bass, paste0(base, ".rds"))

  p_clusters <- plot_spatial(df_bass, paste0("BASS domains (", dataset_tag, ", ", n_pcs, " PCs)"))
  ggsave(paste0(base, "_clusters.png"), p_clusters, width = 6, height = 5, dpi = 300)

  if (!is.null(z_vec)) {
    df_true <- tibble(x = xy_mat[,1], y = xy_mat[,2], z = z_vec)
    p_true <- plot_spatial(df_true, paste0("True labels (", dataset_tag, ")"))
    ggsave(paste0(base, "_true.png"), p_true, width = 6, height = 5, dpi = 300)
  }

  list(df_bass = df_bass, ari = ari_val, bass = bass)
}

# ---- Execute: 10 PCs -----------------------------------------------------------
res10 <- run_bass_with_pcs(
  Y        = Y,
  n_pcs    = 10,
  count_t  = count_t,
  xy_mat   = xy_mat,
  z_vec    = z,
  dataset_tag = "Lung"
)

cat("\n========== SUMMARY (ARI excludes NA) ==========\n")
cat(sprintf("10 PCs : ARI = %.3f\n", ifelse(is.na(res10$ari), NaN, res10$ari)))

end.time_file <- Sys.time()
file_time <- end.time_file - start.time_file
print(file_time)
# ================================================================================