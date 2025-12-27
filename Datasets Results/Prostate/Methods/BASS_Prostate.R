# ====================== BASS on Prostate: 3PC & 10PC =========================
# Counts first, then inject top PCs; ARI computed only where z is not NA
# -----------------------------------------------------------------------------
set.seed(42)

# ---- Libraries --------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(mclust)
library(BASS)
library(Seurat)
library(readr)

# ---- Paths ------------------------------------------------------------------
setwd("~/Desktop/DP-RST-workload/Datasets Results")
out_dir <- "Prostate/Results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load prostate data -----------------------------------------------------
data_dir_prostate_count <- "Prostate/Space_Ranger_Data_prostate/filtered_feature_bc_matrix"
expression_matrix <- Read10X(data.dir = data_dir_prostate_count)  # genes x spots

tissue_positions_list <- readr::read_csv(
  "Prostate/Space_Ranger_Data_prostate/spatial/tissue_positions_list.csv",
  col_names = FALSE
)

# iIMPACT objects: Y (spots x PCs), z (true labels, may contain NA), loc (spots x {x,y}), G
load("Prostate/10x_prostate_cancer_iIMPACT_data.RData")

# ---- Align & build matrices -------------------------------------------------
rownames(loc) <- colnames(expression_matrix)      # ensure loc rows align to barcodes

count_t <- as.matrix(expression_matrix)           # genes x spots
xy_mat  <- as.matrix(loc[, c("x", "y")])          # spots x 2

stopifnot(ncol(count_t) == nrow(xy_mat))
stopifnot(identical(colnames(count_t), rownames(xy_mat)))

# ---- Helper: plotting -------------------------------------------------------
plot_spatial <- function(df, title_text) {
  ggplot(df, aes(x = x, y = y, color = factor(z))) +
    geom_point(size = 0.4) +
    coord_equal() +
    labs(title = title_text, color = "Domain") +
    theme_minimal()
}

# ---- Runner (NA-safe ARI) ---------------------------------------------------
run_bass_with_pcs <- function(Y, n_pcs, count_t, xy_mat, z_vec = NULL, dataset_tag = "Prostate") {
  message("Running BASS with ", n_pcs, " PCs on ", dataset_tag, " ...")
  
  bass <- createBASSObject(
    X  = list(count_t),        # features x cells (genes x spots)
    xy = list(xy_mat),         # cells/spots x 2
    C  = 10, R  = 3,
    beta_method = "SW"
  )
  
  stopifnot(ncol(Y) >= n_pcs)
  pc_mat <- t(as.matrix(Y[, 1:n_pcs, drop = FALSE]))   # n_pcs x spots
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
  
  # ---- ARI excluding NA in true labels -------------------------------------
  ari_val <- NA_real_
  if (!is.null(z_vec)) {
    stopifnot(length(z_vec) == length(zlabels))
    mask <- !is.na(z_vec)
    if (any(mask)) {
      truth_int <- as.integer(as.factor(z_vec[mask]))
      ari_val <- adjustedRandIndex(zlabels[mask], truth_int)
      message("ARI (excluding NA) = ", ari_val)
    } else {
      message("No non-NA labels in z; ARI not computed.")
    }
  }
  
  # Save + plot
  rds_path <- file.path(out_dir, paste0("BASS_Prostate_", n_pcs, "PCs.rds"))
  saveRDS(df_bass, rds_path)
  print(plot_spatial(df_bass, paste0("BASS domains (Prostate, ", n_pcs, " PCs)")))
  
  list(df_bass = df_bass, ari = ari_val, rds_path = rds_path, bass = bass)
}

# ---- Runs: 3 PCs and 10 PCs -------------------------------------------------
res3  <- run_bass_with_pcs(Y, n_pcs = 3,  count_t, xy_mat, z, "Prostate")
res10 <- run_bass_with_pcs(Y, n_pcs = 10, count_t, xy_mat, z, "Prostate")

# ---- Console summary --------------------------------------------------------
cat("\n========== SUMMARY (ARI excludes NA) ==========\n")
cat(sprintf("3 PCs  : ARI = %s\n",  ifelse(is.na(res3$ari), "NA", sprintf("%.6f", res3$ari))))
cat(sprintf("10 PCs : ARI = %s\n",  ifelse(is.na(res10$ari), "NA", sprintf("%.6f", res10$ari))))
cat("Saved RDS:\n")
cat("  - ", res3$rds_path,  "\n", sep = "")
cat("  - ", res10$rds_path, "\n", sep = "")