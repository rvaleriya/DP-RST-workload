# ======================== BASS on STARmap SCE: 3PC & 10PC ====================
# -----------------------------------------------------------------------------
set.seed(42)

# ---- Libraries --------------------------------------------------------------
library(SingleCellExperiment)
library(BASS)
library(ggplot2)
library(dplyr)
library(mclust)

# ---- Paths ------------------------------------------------------------------
setwd("~/Desktop/DP-RST-workload/Datasets Results/Starmap_Mouse")
sce_path <- "starmap_sce_with_pca.rds"
out_dir  <- file.path(getwd(), "Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load SCE ---------------------------------------------------------------
sce <- readRDS(sce_path)

# Basic checks (simple and explicit)
stopifnot("PCA" %in% reducedDimNames(sce))
stopifnot(all(c("x","y") %in% colnames(colData(sce))))
stopifnot("ground_truth" %in% colnames(colData(sce)))

# Counts (genes x cells) and coordinates (cells x 2)
count_mat <- as.matrix(counts(sce))              # genes x cells
coords    <- cbind(x = colData(sce)$x,
                   y = colData(sce)$y)           # cells x 2
barcodes  <- colnames(count_mat)
rownames(coords) <- barcodes

# ---- Helper: plotting -------------------------------------------------------
plot_spatial <- function(df, color_col, title_text) {
  df$.grp <- factor(df[[color_col]])
  ggplot(df, aes(x = x, y = y, color = .grp)) +
    geom_point(size = 0.6, alpha = 0.9) +
    coord_fixed() +
    labs(title = title_text, color = color_col, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

# ---- Runner -----------------------------------------------------------------
run_bass_once_sce <- function(sce, n_pcs) {
  message(sprintf(">> Running BASS with top %d PCs (counts for base X) ...", n_pcs))
  
  # Create BASS object with counts (features x cells)
  bass <- createBASSObject(
    X = list(as.matrix(counts(sce))),     # genes x cells
    xy = list(coords),                    # cells x 2; rownames already set
    C = 20, R = 7,
    beta_method = "SW"
  )
  
  # Inject top n_pcs PCs from SCE PCA (cells x n_pcs) -> transpose to features x cells
  pcs_cells <- reducedDim(sce, "PCA")[, seq_len(n_pcs), drop = FALSE]  # cells x n_pcs
  pcs_feat  <- t(pcs_cells)                                            # n_pcs x cells
  colnames(pcs_feat) <- barcodes                                       # align explicitly
  bass@X_run <- pcs_feat
  
  # Run + postprocess
  bass <- BASS.run(bass)
  bass <- BASS.postprocess(bass, adjustLS = TRUE)
  
  # Labels and evaluation
  zlabels <- as.integer(bass@results$z[[1]])
  truth   <- as.integer(as.factor(colData(sce)$ground_truth))
  stopifnot(length(zlabels) == length(truth))
  ari <- adjustedRandIndex(truth, zlabels)
  
  # Output table
  out_tbl <- tibble(
    ID    = barcodes,
    x     = coords[, "x"],
    y     = coords[, "y"],
    label = zlabels,
    truth = colData(sce)$ground_truth
  )
  
  # Plots
  p_pred  <- plot_spatial(out_tbl, "label",
                          sprintf("BASS — Predicted Domains (%d PCs)", n_pcs))
  p_truth <- plot_spatial(out_tbl, "truth", "Ground Truth — Reference View")
  print(p_pred); print(p_truth)
  
  # Save
  rds_path <- file.path(out_dir, sprintf("BASS_STARmap_%dPCs.rds", n_pcs))
  saveRDS(out_tbl, rds_path)
  
  list(labels_tbl = out_tbl, ari = ari, rds_path = rds_path, bass = bass)
}

# ---- Runs: 3 PCs and 10 PCs -------------------------------------------------
res_3pc  <- run_bass_once_sce(sce, n_pcs = 3)
res_10pc <- run_bass_once_sce(sce, n_pcs = 10)

# ---- Console summary --------------------------------------------------------
cat("\n========== SUMMARY ==========\n")
cat(sprintf("3 PCs  : ARI = %.6f\n",  res_3pc$ari))
cat(sprintf("10 PCs : ARI = %.6f\n",  res_10pc$ari))
cat("Saved RDS:\n")
cat("  - ", res_3pc$rds_path,  "\n", sep = "")
cat("  - ", res_10pc$rds_path, "\n", sep = "")