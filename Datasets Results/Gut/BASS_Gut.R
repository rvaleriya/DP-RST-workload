# ====================== BASS on Intestine (Swiss roll): 3PC & 10PC ======================
# ----------------------------------------------------------------------------------------
set.seed(42)

# ---- Libraries ------------------------------------------------------------------------
library(BASS)
library(Seurat)
library(dplyr)
library(ggplot2)
library(mclust)
library(readr)
library(scCustomize)

# ---- Paths ----------------------------------------------------------------------------
out_dir <- "Gut/Results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load precomputed PCs + labels ----------------------------------------------------
gut_df_wt_muscle <- readRDS("Gut/gut_df_wt_muscle.rds")
# Expect: rownames = barcodes; columns PC1.., x, y, z_man, z

# ---- Load counts + coords as you already prepared ------------------------------------
### If you haven't created the filtered Seurat object in this session, run your block first:
data_dir <- 'Gut/Space_Ranger_Data_Gut'
expression_matrix <- Read10X_h5_GEO(data_dir = data_dir)
load("Gut/swiss_roll_wt_muscle_finaltouches1.RData")
loc = swiss_roll_wt_muscle_finaltouches1[, 4:5]
colnames(loc) <- c("row", "col")
seurat_object = CreateSeuratObject(counts = expression_matrix, meta.data = loc)
meta <- seurat_object@meta.data
filtered_meta <- meta[!is.na(meta$row) & !is.na(meta$col), ]
seurat_object_filtered <- subset(seurat_object, cells = rownames(filtered_meta))

# For this script we assume you already have:
stopifnot(exists("seurat_object_filtered"))

# ---- Build matrices for BASS ----------------------------------------------------------
# Counts: genes x spots
count_mat <- as.matrix(GetAssayData(seurat_object_filtered, assay = "RNA", slot = "counts"))
barcodes  <- colnames(count_mat)

# Coordinates from meta.data (rename to x,y)
meta_f <- seurat_object_filtered@meta.data
stopifnot(all(c("row", "col") %in% colnames(meta_f)))
xy_mat <- as.matrix(meta_f[, c("row", "col")])
colnames(xy_mat) <- c("x", "y")
# Ensure ordering matches barcodes
xy_mat <- xy_mat[barcodes, , drop = FALSE]
stopifnot(identical(rownames(xy_mat), barcodes))

# ---- Align PCs and labels from gut_df_wt_muscle --------------------------------------
# Keep only barcodes present in the counts
gut_sub <- gut_df_wt_muscle[barcodes, , drop = FALSE]

# Extract a simple PC matrix accessor
get_pc_mat <- function(df, n_pcs) {
  pc_cols <- paste0("PC", seq_len(n_pcs))
  stopifnot(all(pc_cols %in% colnames(df)))
  as.matrix(df[, pc_cols, drop = FALSE])  # spots x n_pcs
}

# True labels vector (may contain NA)
z_vec <- if ("z" %in% colnames(gut_sub)) gut_sub[["z"]] else NULL

# ---- Helper: plotting ---------------------------------------------------------------
plot_spatial <- function(df, title_text) {
  ggplot(df, aes(x = x, y = y, color = factor(z))) +
    geom_point(size = 0.4) +
    coord_equal() +
    labs(title = title_text, color = "Domain") +
    theme_minimal()
}

# ---- Runner (counts -> BASS; inject PCs; ARI excludes NA) ---------------------------
run_bass_intestine <- function(count_mat, xy_mat, pc_source_df, n_pcs, z_vec = NULL) {
  message("Running BASS with ", n_pcs, " PCs on Intestine ...")
  
  # Create BASS object from counts (features x cells)
  bass <- createBASSObject(
    X  = list(count_mat),    # genes x spots
    xy = list(xy_mat),       # spots x 2 (x,y), rownames = barcodes
    C  = 15, R  = 5,
    beta_method = "SW"
  )
  
  # Inject top n PCs (spots x n_pcs) -> transpose to n_pcs x spots
  pcs_spots <- get_pc_mat(pc_source_df, n_pcs)      # spots x n_pcs
  pcs_feat  <- t(pcs_spots)                         # n_pcs x spots
  colnames(pcs_feat) <- colnames(count_mat)         # align columns
  bass@X_run <- pcs_feat
  
  # Run + postprocess
  bass <- BASS.run(bass)
  bass <- BASS.postprocess(bass, adjustLS = TRUE)
  
  # Labels
  z_pred <- as.integer(bass@results$z[[1]])
  df_out <- tibble(
    barcode = colnames(count_mat),
    x = xy_mat[, "x"],
    y = xy_mat[, "y"],
    z = z_pred
  )
  
  # ARI (exclude NA in truth)
  ari_val <- NA_real_
  if (!is.null(z_vec)) {
    stopifnot(length(z_vec) == length(z_pred))
    mask <- !is.na(z_vec)
    if (any(mask)) {
      truth_int <- as.integer(as.factor(z_vec[mask]))
      ari_val <- adjustedRandIndex(z_pred[mask], truth_int)
      message("ARI (excluding NA) = ", ari_val)
    } else {
      message("No non-NA labels in z; ARI not computed.")
    }
  }
  
  list(df = df_out, ari = ari_val, bass = bass)
}

# ---- Runs: 3 PCs and 10 PCs --------------------------------------------------------
res3  <- run_bass_intestine(count_mat, xy_mat, gut_sub, n_pcs = 3,  z_vec = z_vec)
res10 <- run_bass_intestine(count_mat, xy_mat, gut_sub, n_pcs = 10, z_vec = z_vec)

# ---- Save + Plot --------------------------------------------------------------------
rds3  <- file.path(out_dir, "BASS_Intestine_3PCs.rds")
rds10 <- file.path(out_dir, "BASS_Intestine_10PCs.rds")
saveRDS(res3$df,  rds3)
saveRDS(res10$df, rds10)

print(plot_spatial(res3$df,  "BASS domains (Intestine, 3 PCs)"))
print(plot_spatial(res10$df, "BASS domains (Intestine, 10 PCs)"))

cat("\n========== SUMMARY (Intestine; ARI excludes NA) ==========\n")
cat(sprintf("3 PCs  : ARI = %s\n",  ifelse(is.na(res3$ari),  "NA", sprintf("%.6f", res3$ari))))
cat(sprintf("10 PCs : ARI = %s\n",  ifelse(is.na(res10$ari), "NA", sprintf("%.6f", res10$ari))))
cat("Saved RDS:\n")
cat("  - ", rds3,  "\n", sep = "")
cat("  - ", rds10, "\n", sep = "")