# ====================== k-means on osmFISH: 3PC & 10PC =======================

set.seed(42)

# ---- Libraries --------------------------------------------------------------
library(readr)
library(dplyr)
library(ggplot2)
library(mclust)   

# ---- Paths ------------------------------------------------------------------
setwd("~/Desktop/DP-RST-workload/Datasets Results/osmFISH")
in_csv   <- "osmfish_pcs_coords_labels.csv"
out_dir  <- file.path(getwd(), "Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Data -------------------------------------------------------------------
osmfish <- read_csv(in_csv, show_col_types = FALSE)

# ---- Helpers ----------------------------------------------------------------
plot_spatial <- function(df, color_col, title_text) {
  stopifnot(color_col %in% names(df))
  df$.grp <- factor(df[[color_col]])
  
  ggplot(df, aes(x = x, y = y, color = .grp)) +
    geom_point(size = 0.6, alpha = 0.9) +
    coord_fixed() +
    labs(title = title_text, color = color_col, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

run_kmeans_once <- function(df, n_pcs) {
  message(sprintf(">> Running k-means (K=11) with %d PCs ...", n_pcs))
  
  pc_cols <- paste0("PC", seq_len(n_pcs))
  stopifnot(all(pc_cols %in% colnames(df)))
  X <- as.matrix(df[, pc_cols])
  
  # k-means clustering
  km <- kmeans(X, centers = 11)
  labels <- as.integer(km$cluster)
  
  # ARI vs ground truth
  ref <- as.integer(as.factor(df$ground_truth))
  stopifnot(length(ref) == length(labels))
  ari <- adjustedRandIndex(ref, labels)
  
  out_tbl <- tibble(
    ID     = df$ID,
    x      = df$x,
    y      = df$y,
    label  = labels,
    truth  = df$ground_truth
  )
  
  # Plot directly in session
  p_pred  <- plot_spatial(out_tbl, "label",
                          sprintf("k-means (K=11), %d PCs — Predicted Clusters", n_pcs))
  p_truth <- plot_spatial(out_tbl, "truth",
                          "Ground Truth — Reference View")
  
  print(p_pred)
  print(p_truth)
  
  # Save labels to RDS
  rds_path <- file.path(out_dir, sprintf("KMeans_osmFISH_%dPCs.rds", n_pcs))
  saveRDS(out_tbl, rds_path)
  
  list(labels_tbl = out_tbl, ari = ari, rds_path = rds_path)
}

# ---- Runs: 3 PCs and 10 PCs -------------------------------------------------
res_3pc  <- run_kmeans_once(osmfish, n_pcs = 3)
res_10pc <- run_kmeans_once(osmfish, n_pcs = 10)

# ---- Console summary --------------------------------------------------------
cat("\n========== SUMMARY ==========\n")
cat(sprintf("3 PCs  : ARI = %.6f\n",  res_3pc$ari))
cat(sprintf("10 PCs : ARI = %.6f\n",  res_10pc$ari))
cat("Saved RDS:\n")
cat("  - ", res_3pc$rds_path,  "\n", sep = "")
cat("  - ", res_10pc$rds_path, "\n", sep = "")