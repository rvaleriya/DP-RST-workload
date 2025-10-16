# ---------------------------------- Setup -------------------------------------
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("zhengli09/BASS")

library(ggplot2)
library(dplyr)
library(purrr)
library(mclust)
library(tidyr)
library(stringr)
library(BASS)

set.seed(42)
setwd("~/Desktop/DP-RST-workload/Datasets Results")

# Data objects expected from the .RData file:
#   count (cells x genes), loc (spots x 2), Y (spots x PCs), ground_truth_df (x,y,type)
load("Brain/DFPLC_151510_data_for_iIMPACT.RData")

# --------------------------------- BASS run -----------------------------------
# Align coordinates to barcodes used in expression matrix
#    After transpose, t(count) is genes x cells; its colnames are spot barcodes.
rownames(loc) <- colnames(t(count))

# Transpose counts: genes x cells
count_t <- t(as.matrix(count))
xy_mat  <- as.matrix(loc[, c("x", "y")])

# Helper function to run BASS on selected PCs
run_bass_with_pcs <- function(Y, n_pcs, count_t, xy_mat, ground_truth_df, out_prefix) {
  message("Running BASS with ", n_pcs, " PCs ...")
  
  # Create BASS object
  bass <- createBASSObject(
    X = list(count_t),
    xy = list(xy_mat),
    C = 20, R = 7,
    beta_method = "SW"
  )
  
  # Inject top n_pcs PCs
  pc_mat <- t(as.matrix(Y[, 1:n_pcs, drop = FALSE]))
  bass@X_run <- pc_mat
  
  # Run + postprocess
  bass <- BASS.run(bass)
  bass <- BASS.postprocess(bass, adjustLS = TRUE)
  
  # Extract domain labels
  zlabels <- bass@results$z[[1]]
  
  # Create dataframe with coords and domains
  df_bass <- tibble(
    barcode = colnames(count_t),
    x = xy_mat[, 1],
    y = xy_mat[, 2],
    z = zlabels
  )
  
  # Merge with ground truth
  df_eval <- df_bass %>%
    inner_join(ground_truth_df, by = c("x", "y")) %>%
    filter(!is.na(type))
  
  # Compute ARI
  ari <- adjustedRandIndex(df_eval$z, df_eval$type)
  message("ARI with ", n_pcs, " PCs: ", ari)
  
  # Save outputs
  labels_out <- df_bass %>%
    select(barcode, x, y, domain = z)
  
  saveRDS(labels_out, paste0("Brain/Results/BASS_Brain_", n_pcs, "PCs.rds"))
  # write.csv(labels_out, paste0("Brain/Results/BASS_Brain_", n_pcs, "PCs.csv"), row.names = FALSE)
  # writeLines(as.character(ari), paste0("Brain/Results/BASS_Brain_ARI_", n_pcs, "PCs.txt"))
  
  # Return for plotting if desired
  list(bass = bass, df_bass = df_bass, ari = ari)
}

# ---------------- Run with 3 PCs ----------------
res3 <- run_bass_with_pcs(Y, n_pcs = 3, count_t, xy_mat, ground_truth_df, "Brain")

# Plot
ggplot(res3$df_bass, aes(x, y, color = factor(z))) +
  geom_point(size = 0.4) +
  coord_equal() +
  labs(title = "BASS domains (3 PCs)", color = "Domain") +
  theme_minimal()

# ---------------- Run with 10 PCs ----------------
res10 <- run_bass_with_pcs(Y, n_pcs = 10, count_t, xy_mat, ground_truth_df, "Brain")

# Plot
ggplot(res10$df_bass, aes(x, y, color = factor(z))) +
  geom_point(size = 0.4) +
  coord_equal() +
  labs(title = "BASS domains (10 PCs)", color = "Domain") +
  theme_minimal()