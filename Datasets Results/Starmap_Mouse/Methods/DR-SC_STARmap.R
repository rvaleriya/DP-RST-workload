##### DR-SC on STARmap (simple version, fixed coordinate names) ###############
set.seed(42)

library(SingleCellExperiment)
library(Seurat)
library(DR.SC)
library(dplyr)
library(ggplot2)
library(mclust)
library(readr)
library(FNN)

# ---- Paths ------------------------------------------------------------------
setwd("~/Desktop/DP-RST-workload/Datasets Results/Starmap_Mouse")
sce_rds <- "starmap_sce_with_pca.rds"
out_dir <- file.path(getwd(), "Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load SCE and create Seurat object --------------------------------------
sce <- readRDS(sce_rds)

# counts matrix
if ("counts" %in% assayNames(sce)) {
  mat <- assay(sce, "counts")
} else {
  stop("No 'counts' assay found in the SCE object.")
}

# coordinates and metadata
cd <- as.data.frame(colData(sce))

# --- DR.SC expects 'row' and 'col' names, not 'x'/'y' -------------------
if (all(c("x","y") %in% names(cd))) {
  meta <- cd[, c("x","y", intersect("ground_truth", names(cd))), drop = FALSE]
  names(meta)[1:2] <- c("row","col")
} else if (all(c("row","col") %in% names(cd))) {
  meta <- cd[, c("row","col", intersect("ground_truth", names(cd))), drop = FALSE]
} else {
  stop("Need spatial coordinates in colData(sce) as 'x,y' or 'row,col'.")
}
rownames(meta) <- colnames(sce)

seu <- CreateSeuratObject(counts = mat, meta.data = meta)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)

# ---- Helper for plotting ----------------------------------------------------
plot_spatial <- function(df, color_col, title_text) {
  df$.grp <- factor(df[[color_col]])
  ggplot(df, aes(x = col, y = row, color = .grp)) +   # note row/col for plotting
    geom_point(size = 0.6, alpha = 0.9) +
    coord_fixed() +
    labs(title = title_text, color = color_col, x = "col", y = "row") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

# ---- Build adjacency manually (radius-based) --------------------------------
pos <- as.matrix(seu@meta.data[, c("row","col")])

# Heuristic radius cap: median distance to the 6th neighbor × 1.2
# (ensures ~ ≥6 neighbors for most points on irregular STARmap layouts)
nn6 <- get.knn(pos, k = 6)$nn.dist[, 6]
rad_cap <- median(nn6, na.rm = TRUE) * 1.2

# Try to target median neighbors ~ 6–8
Adj_mat <- DR.SC::getAdj_auto(pos, lower.med = 6, upper.med = 8, radius.upper = rad_cap)

# ---- Prepare expression X (log-normalized var.features) ------------------
if (inherits(seu@assays$RNA, "Assay5")) {
  var.features <- seu@assays$RNA@meta.data$var.features
  var.features <- var.features[!is.na(var.features)]
  dat <- Seurat::GetAssayData(seu, assay = "RNA", layer = "data")
  X <- Matrix::t(dat[var.features, ])
} else {
  var.features <- seu@assays$RNA@var.features
  X <- Matrix::t(seu[["RNA"]]@data[var.features, ])
}

# Simple DR-SC_fit wrapper (expects: X = cells×features; rownames(X)=cell IDs;
# Adj_mat = symmetric sparse adjacency with dim = nrow(X); plot_spatial() defined)
run_drsc_once <- function(seu_obj, q_dim, K = 7) {
  message(sprintf(">> DR.SC_fit | q=%d, K=%d", q_dim, K))
  
  fit <- DR.SC::DR.SC_fit(
    X      = X,
    K      = K,
    q      = q_dim,
    Adj_sp = Adj_mat,
    verbose = FALSE
  )
  
  # Cluster labels (one per row of X)
  lab <- as.integer(fit$Objdrsc[[1]]$cluster)
  
  # Cell IDs to align everything
  ids <- if (!is.null(rownames(X))) rownames(X) else colnames(seu_obj)
  
  # Build output table aligned by IDs
  md <- seu_obj@meta.data
  md$ID <- rownames(md)
  
  out_tbl <- tibble::tibble(
    ID   = ids,
    row  = md[ids, "row"],
    col  = md[ids, "col"],
    label = lab
  )
  
  # Optional ARI if ground truth exists
  if ("ground_truth" %in% colnames(md)) {
    truth <- md[ids, "ground_truth"]
    out_tbl$truth <- truth
    keep <- !is.na(out_tbl$label) & !is.na(truth)
    ari  <- if (any(keep)) mclust::adjustedRandIndex(
      as.integer(as.factor(truth[keep])), out_tbl$label[keep]
    ) else NA_real_
  } else {
    ari <- NA_real_
  }
  
  # Plots
  p_pred <- plot_spatial(out_tbl, "label",
                         sprintf("DR-SC (K=%d, q=%d) — Predicted Clusters", K, q_dim))
  print(p_pred)
  if ("truth" %in% names(out_tbl)) {
    p_truth <- plot_spatial(out_tbl, "truth", "Ground Truth — Reference View")
    print(p_truth)
  }
  
  # Save
  rds_path <- file.path(out_dir, sprintf("DR-SC_STARmap_%dPCs.rds", q_dim))
  saveRDS(out_tbl, rds_path)
  
  list(labels_tbl = out_tbl, ari = ari, rds_path = rds_path)
}

# Run both settings
res_q3  <- run_drsc_once(seu, q_dim = 3)
res_q10 <- run_drsc_once(seu, q_dim = 10)

cat("\n========== SUMMARY ==========\n")
cat(sprintf("q=3  : ARI = %s\n",  ifelse(is.na(res_q3$ari),  "NA", format(res_q3$ari, digits=6))))
cat(sprintf("q=10 : ARI = %s\n",  ifelse(is.na(res_q10$ari), "NA", format(res_q10$ari, digits=6))))
cat("Saved RDS:\n")
cat("  - ", res_q3$rds_path,  "\n", sep = "")
cat("  - ", res_q10$rds_path, "\n", sep = "")