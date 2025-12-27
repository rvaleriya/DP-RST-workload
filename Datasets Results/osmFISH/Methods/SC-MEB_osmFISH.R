# ====================== SC.MEB on osmFISH: 3PC & 10PC =========================

# -----------------------------------------------------------------------------
set.seed(42)

# ---- Libraries --------------------------------------------------------------
library(readr)
library(dplyr)
library(ggplot2)
library(mclust)
library(SC.MEB)
library(FNN)

# ---- Paths ------------------------------------------------------------------
setwd("~/Desktop/DP-RST-workload/Datasets Results/osmFISH")
in_csv   <- "osmfish_pcs_coords_labels.csv"
out_dir  <- file.path(getwd(), "Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Data -------------------------------------------------------------------
starmap <- read_csv(in_csv, show_col_types = FALSE)

# ---- SC.MEB hyper-parameters ------------------------------------------------
beta_grid     <- seq(0, 5, 0.2)
K_set         <- c(11)
parallel_run  <- TRUE
num_core      <- 3
PX            <- TRUE
maxIter_ICM   <- 10
maxIter       <- 50

# ---- Spatial neighborhood ---------------------------------------------------
coords_mat <- as.matrix(starmap[, c("x","y")])
Adj_sp <- getneighborhood_fast(coords_mat, cutoff = 400)

# ---- Helpers ----------------------------------------------------------------
extract_labels <- function(scmeb_fit_singleK) {
  x_obj <- scmeb_fit_singleK$x
  if (is.matrix(x_obj)) as.integer(x_obj[, 1]) else as.integer(x_obj)
}

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

run_scmeb_once <- function(df, n_pcs) {
  message(sprintf(">> Running SC.MEB with %d PCs ...", n_pcs))
  
  pc_cols <- paste0("PC", seq_len(n_pcs))
  stopifnot(all(pc_cols %in% colnames(df)))
  X <- as.matrix(df[, pc_cols])
  
  fit <- SC.MEB(
    y            = X,
    Adj_sp       = Adj_sp,
    beta_grid    = beta_grid,
    K_set        = K_set,
    parallel     = parallel_run,
    num_core     = num_core,
    PX           = PX,
    maxIter_ICM  = maxIter_ICM,
    maxIter      = maxIter
  )
  
  resK   <- fit[, 1]
  labels <- extract_labels(resK)
  
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
                          sprintf("SC.MEB (K=8), %d PCs — Predicted Clusters", n_pcs))
  p_truth <- plot_spatial(out_tbl, "truth",
                          "Ground Truth — Reference View")
  
  print(p_pred)
  print(p_truth)
  
  # Save labels to RDS
  rds_path <- file.path(out_dir, sprintf("SC-MEB_osmFISH_%dPCs.rds", n_pcs))
  saveRDS(out_tbl, rds_path)
  
  list(labels_tbl = out_tbl, ari = ari, rds_path = rds_path)
}

# ---- Runs: 3 PCs and 10 PCs -------------------------------------------------
res_3pc  <- run_scmeb_once(starmap, n_pcs = 3)
res_10pc <- run_scmeb_once(starmap, n_pcs = 10)

# ---- Console summary --------------------------------------------------------
cat("\n========== SUMMARY ==========\n")
cat(sprintf("3 PCs  : ARI = %.6f\n",  res_3pc$ari))
cat(sprintf("10 PCs : ARI = %.6f\n",  res_10pc$ari))
cat("Saved RDS:\n")
cat("  - ", res_3pc$rds_path,  "\n", sep = "")
cat("  - ", res_10pc$rds_path, "\n", sep = "")
