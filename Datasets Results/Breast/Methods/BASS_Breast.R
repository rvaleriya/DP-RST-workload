# ---------------------------------- Setup -------------------------------------
library(ggplot2)
library(dplyr)
library(purrr)
library(mclust)
library(tidyr)
library(stringr)
library(BASS)
library(Seurat)
library(readr)

set.seed(42)
setwd("~/Desktop/DP-RST-workload/Datasets Results")

# ---------------------------------- Load data ---------------------------------
expression_matrix <- Read10X("Breast/Space_Ranger_Data_Breast/filtered_feature_bc_matrix")

tissue_positions_list <- readr::read_csv(
  "Breast/Space_Ranger_Data_Breast/spatial/tissue_positions_list.csv",
  col_names = FALSE
)
# Load iIMPACT objects (Y, z, loc, G)
load("Breast/10x_breast_cancer_iIMPACT_data.RData")

# -------------------------------- Align data ----------------------------------
# Make sure loc rows match expression_matrix columns
rownames(loc) <- colnames(expression_matrix)

# Build matrices for BASS
count_t <- as.matrix(expression_matrix)  # genes x spots (already correct orientation for BASS)
xy_mat  <- as.matrix(loc[, c("x", "y")])

# --------------------------------- BASS run -----------------------------------
run_bass_with_pcs <- function(Y, n_pcs, count_t, xy_mat, z_vec = NULL, dataset_tag = "Breast") {
  message("Running BASS with ", n_pcs, " PCs on ", dataset_tag, " ...")
  
  bass <- createBASSObject(
    X  = list(count_t),
    xy = list(xy_mat),
    C  = 15, R  = 5,
    beta_method = "SW"
  )
  
  pc_mat <- t(as.matrix(Y[, 1:n_pcs, drop = FALSE]))
  bass@X_run <- pc_mat
  
  bass <- BASS.run(bass)
  bass <- BASS.postprocess(bass, adjustLS = TRUE)
  
  zlabels <- bass@results$z[[1]]
  
  df_bass <- tibble(
    barcode = colnames(count_t),
    x = xy_mat[, 1],
    y = xy_mat[, 2],
    z = zlabels
  )
  
  # ARI (if z available)
  ari_val <- NA_real_
  if (!is.null(z_vec)) {
    df_eval <- df_bass %>%
      mutate(type = z_vec) %>%
      filter(!is.na(type))
    ari_val <- adjustedRandIndex(df_eval$z, df_eval$type)
    message("ARI = ", ari_val)
  }
  
  # Save
  dir.create("Breast/Results", showWarnings = FALSE, recursive = TRUE)
  saveRDS(df_bass, paste0("Breast/Results/BASS_Breast_", n_pcs, "PCs.rds"))
  
  list(df_bass = df_bass, ari = ari_val)
}

# ---------------- Run with 3 PCs ----------------
res3 <- run_bass_with_pcs(Y, n_pcs = 3, count_t, xy_mat, z, "Breast")

ggplot(res3$df_bass, aes(x, y, color = factor(z))) +
  geom_point(size = 0.4) +
  coord_equal() +
  labs(title = "BASS domains (Breast, 3 PCs)", color = "Domain") +
  theme_minimal()

# ---------------- Run with 10 PCs ----------------
res10 <- run_bass_with_pcs(Y, n_pcs = 10, count_t, xy_mat, z, "Breast")

ggplot(res10$df_bass, aes(x, y, color = factor(z))) +
  geom_point(size = 0.4) +
  coord_equal() +
  labs(title = "BASS domains (Breast, 10 PCs)", color = "Domain") +
  theme_minimal()