#!/usr/bin/env Rscript
# ==============================================================================
# Lung Xenium Data Processing Pipeline (VUILD96MF)
# ------------------------------------------------------------------------------
# Steps:
#   1. Load Xenium counts + cells + annotations
#   2. Align counts to cells; compute total_counts
#   3. Apply quality control filters:
#        - Keep cells inside tissue boundary
#        - Identify and remove high-count outliers using the Kneedle method
#        - Keep genes detected in >= 3 cells
#   4. Seurat preprocessing:
#        - LogNormalize
#        - Select 2,000 HVGs
#        - Run PCA (15 components)
#   5. Save outputs in multiple formats:
#        - Seurat object (.rds)
#        - Processed RData bundle
#        - AnnData (.h5ad) for Python/Scanpy
#
# Outputs:
#   - Data_VUILD96MF/Lung_xenium_seurat_object_processed.rds
#   - Data_VUILD96MF/Lung_xenium_data_processed.RData
#   - Data_VUILD96MF/Lung_xenium_processed.h5ad
# ==============================================================================

library(Matrix)
library(tidyverse)
library(data.table)
library(sparseMatrixStats)
library(Seurat)
library(sf)
library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)

set.seed(42)

# ---- Paths -------------------------------------------------------------------
setwd("~/Desktop/DP-RST-workload/Datasets Results/Lung_Xenium")

transcript_data <- "Data_VUILD96MF/cell_feature_matrix.h5"
cells_data      <- "Data_VUILD96MF/cells.csv.gz"
annotation_data <- "Data_VUILD96MF/cells_partitioned_by_annotation.csv"
boundary_csv    <- "boundary_outer.csv"

rds_out         <- "Data_VUILD96MF/Lung_xenium_seurat_object_processed.rds"
rdata_out       <- "Data_VUILD96MF/Lung_xenium_data_processed.RData"
h5ad_out        <- "Data_VUILD96MF/Lung_xenium_processed.h5ad"

# ---- Load Xenium counts, cells, annotations ---------------------------------
cat("\n=== LOADING XENIUM DATA ===\n")

expr_list  <- Read10X_h5(transcript_data)          # list, e.g. "Gene Expression"
counts_all <- expr_list[[1]]
if (!inherits(counts_all, "dgCMatrix")) counts_all <- as(counts_all, "dgCMatrix")
cat(sprintf("Counts: %d genes × %d cells\n", nrow(counts_all), ncol(counts_all)))

cells <- fread(cells_data) |> as_tibble()
annotations <- read_csv(annotation_data, show_col_types = FALSE) |>
  filter(str_detect(sample, "VUILD96MF"))

# De-duplicate annotations by cell_id (keep first)
dup_ids <- names(which(table(annotations$cell_id) > 1))
annotations_unique <- annotations |>
  filter(!(cell_id %in% dup_ids)) |>
  select(cell_id, Annotation_Type)

# Left join annotations → cells_annotated (keep format)
cells_annotated <- cells %>%
  left_join(
    annotations_unique %>% select(cell_id, Annotation_Type),
    by = "cell_id"
  )

cat(sprintf("Cells after join: %d\n", nrow(cells_annotated)))

# ---- Align counts to cell metadata; compute total_counts ---------------------
cat("\n=== ALIGNING COUNTS TO CELLS & COMPUTING TOTAL COUNTS ===\n")
inter <- intersect(colnames(counts_all), cells_annotated$cell_id)

cells_annotated <- cells_annotated %>% filter(cell_id %in% inter)
cells_annotated <- cells_annotated[match(inter, cells_annotated$cell_id), ]
counts_all      <- counts_all[, inter, drop = FALSE]

cells_annotated$total_counts <- as.numeric(Matrix::colSums(counts_all))
cat(sprintf("Aligned: %d cells | median total_counts: %.0f\n",
            ncol(counts_all), median(cells_annotated$total_counts)))

# ---- Read boundary and build polygons ---------------------------------------
cat("\n=== READING BOUNDARY & BUILDING POLYGONS ===\n")
bnd <- read.csv(boundary_csv, check.names = FALSE)

# Ensure the polygon is closed (first point equals last point)
bmat <- as.matrix(bnd[, c("x","y")])
if (!all(bmat[1, ] == bmat[nrow(bmat), ])) {
  bmat <- rbind(bmat, bmat[1, ])
}

# Build a single polygon (give it an internal id=1 for convenience)
poly_single <- sf::st_polygon(list(bmat))
polys_sf    <- sf::st_sf(polygon_id = 1L, geometry = sf::st_sfc(poly_single))

cat("Constructed 1 boundary polygon from x,y.\n")

# ---- Spatial filter: inside boundary ----------------------------------------
cat("\n=== FILTERING CELLS INSIDE BOUNDARY ===\n")
pts_sf <- sf::st_as_sf(cells_annotated,
                       coords = c("x_centroid","y_centroid"),
                       remove = FALSE, crs = sf::NA_crs_)
inside <- lengths(sf::st_intersects(pts_sf, polys_sf)) > 0
cat(sprintf("Inside boundary: %d / %d cells\n", sum(inside), nrow(cells_annotated)))

# ---- Kneedle upper cutoff on total_counts (baseline keep_upper) -------------
cat("\n=== KNEEDLE UPPER CUTOFF ON TOTAL COUNTS (BASELINE) ===\n")
tc_in <- cells_annotated$total_counts[inside]
tc_sorted <- sort(tc_in, decreasing = TRUE)
n <- length(tc_sorted)

x  <- seq_len(n)
xn <- (x - x[1]) / (x[n] - x[1])
yn <- (tc_sorted - tc_sorted[n]) / (tc_sorted[1] - tc_sorted[n])
d  <- abs(-1 * xn - 1 * yn + 1) / sqrt(2)   # distance to line (0,1)->(1,0)
knee_idx  <- which.max(d)
tc_kneedle <- tc_sorted[knee_idx]

keep_upper <- cells_annotated$total_counts <= tc_kneedle
cat(sprintf("Kneedle cutoff: %.2f | Pass: %d (%.2f%%)\n",
            tc_kneedle, sum(keep_upper), 100*mean(keep_upper)))

# ==== ENHANCED XENIUM QC (PF-AWARE) ===========================================
cat("\n=== XENIUM QC (PF-AWARE) ===\n")

eps <- 1e-9

cells_annotated <- cells_annotated %>%
  mutate(
    total_counts_safe = pmax(total_counts, eps),
    control_ratio     = pmin(pmax(control_probe_counts       / total_counts_safe, 0), 1),
    unassigned_ratio  = pmin(pmax(unassigned_codeword_counts / total_counts_safe, 0), 1),
    assignment_ratio  = pmin(pmax(transcript_counts          / total_counts_safe, 0), 1),
    nucleus_cell_ratio = nucleus_area / pmax(cell_area, eps)
  )

# Thresholds 
min_transcripts      <- 18
min_cell_area        <- 20
min_assignment_ratio <- 0.50

# Area: MAD rule + PF hard cap
log_area        <- log10(cells_annotated$cell_area + 1)
median_log_area <- median(log_area, na.rm = TRUE)
mad_log_area    <- mad(log_area, na.rm = TRUE)
max_area_mad    <- 10^(median_log_area + 2.5 * mad_log_area) - 1
max_cell_area_final <- min(max_area_mad, 1000)  # liberal PF cap (≈ up to 1000 µm^2)

# Chemistry/morphology bounds
max_control_ratio    <- 0.10
max_unassigned_ratio <- 0.12
min_nucleus_ratio    <- 0.03
max_nucleus_ratio    <- 0.85

# Primary per-metric gates 
keep_transcripts <- cells_annotated$transcript_counts >= min_transcripts
keep_area <- (cells_annotated$cell_area >= min_cell_area) &
  (cells_annotated$cell_area <= max_cell_area_final)
keep_control    <- cells_annotated$control_ratio    <= max_control_ratio
keep_unassigned <- cells_annotated$unassigned_ratio <= max_unassigned_ratio
keep_nucleus    <- (cells_annotated$nucleus_cell_ratio >= min_nucleus_ratio) &
  (cells_annotated$nucleus_cell_ratio <= max_nucleus_ratio)
keep_assignment <- cells_annotated$assignment_ratio >= min_assignment_ratio

# Doublets on log–log residuals (97th percentile)
cells_annotated <- cells_annotated %>%
  mutate(
    log_area_pf = log10(cell_area + 1),
    log_tx_pf   = log10(transcript_counts + 1)
  )
fit_pf   <- lm(log_tx_pf ~ log_area_pf, data = cells_annotated)
res_pf   <- as.numeric(scale(residuals(fit_pf)))
dbl_score <- pmax(0, res_pf) + pmax(0, as.numeric(scale(cells_annotated$cell_area)))
dbl_thr   <- quantile(dbl_score, 0.97, na.rm = TRUE)
likely_doublets <- dbl_score > dbl_thr

# Relax Kneedle using morphology 
keep_upper <- keep_upper | (
  (cells_annotated$cell_area <= max_cell_area_final) &
    !likely_doublets &
    (cells_annotated$nucleus_cell_ratio >= min_nucleus_ratio) &
    (cells_annotated$nucleus_cell_ratio <= max_nucleus_ratio)
)

# Combine all 
keep_final <- inside & keep_upper &
  keep_transcripts & keep_area &
  keep_control & keep_unassigned &
  keep_nucleus & keep_assignment &
  !likely_doublets

cat("\n--- Filtering Summary (PF) ---\n")
cat(sprintf("Inside boundary        : %d (%.2f%%)\n", sum(inside), 100*mean(inside)))
cat(sprintf("Count <= Kneedle/relax : %d (%.2f%%)\n", sum(keep_upper), 100*mean(keep_upper)))
cat(sprintf("Min transcripts (>= %d): %d (%.2f%%)\n", min_transcripts, sum(keep_transcripts), 100*mean(keep_transcripts)))
cat(sprintf("Cell area OK           : %d (%.2f%%)\n", sum(keep_area), 100*mean(keep_area)))
cat(sprintf("Control probes OK      : %d (%.2f%%)\n", sum(keep_control), 100*mean(keep_control)))
cat(sprintf("Unassigned OK          : %d (%.2f%%)\n", sum(keep_unassigned), 100*mean(keep_unassigned)))
cat(sprintf("Nucleus ratio OK       : %d (%.2f%%)\n", sum(keep_nucleus), 100*mean(keep_nucleus)))
cat(sprintf("Assignment ratio OK    : %d (%.2f%%)\n", sum(keep_assignment), 100*mean(keep_assignment)))
cat(sprintf("Not likely doublets    : %d (%.2f%%)\n", sum(!likely_doublets), 100*mean(!likely_doublets)))
cat(sprintf("FINAL cells kept       : %d (%.2f%%)\n", sum(keep_final), 100*mean(keep_final)))

# Persist per-cell flags for audit/repro
cells_annotated$keep_final      <- keep_final
cells_annotated$keep_upper_rel  <- keep_upper
cells_annotated$likely_doublet  <- likely_doublets

# Quick visualization
cells_annotated$qc_status <- ifelse(keep_final, "Kept", "Excluded")
ggplot(cells_annotated, aes(x = x_centroid, y = y_centroid)) +
  geom_point(aes(color = qc_status), size = 0.3, alpha = 0.7) +
  scale_color_manual(values = c("Excluded" = "red", "Kept" = "grey70")) +
  coord_fixed() + theme_minimal() + theme(legend.position = "top") +
  labs(title = "Spatial Map of Excluded vs Kept Cells", x = "X", y = "Y")


# ---- Seurat preprocessing ----------------------------------------------------
cat("\n=== SEURAT PREPROCESSING (LogNormalize -> 2,000 HVGs -> 15 PCs) ===\n")

spatial_keep <- cells_annotated[keep_final, , drop = FALSE]
counts_keep  <- counts_all[, spatial_keep$cell_id, drop = FALSE]

so <- CreateSeuratObject(
  counts    = counts_keep,
  meta.data = spatial_keep |> as.data.frame(),
  project   = "Lung_Xenium_VUILD96MF",
  assay     = "RNA",
  min.cells = 0, min.features = 0
)

so <- NormalizeData(so, normalization.method = "LogNormalize",
                    scale.factor = 1e4, verbose = FALSE)
so <- FindVariableFeatures(so, selection.method = "vst",
                           nfeatures = min(2000, nrow(counts_keep)), 
                           verbose = FALSE)
hvg_genes <- VariableFeatures(so)
so <- ScaleData(so, features = hvg_genes, verbose = FALSE)
so <- RunPCA(so, features = hvg_genes, npcs = 15, verbose = FALSE)
cat(sprintf("Computed HVGs: %d | PCs: %d\n", length(hvg_genes), 15))

# ---- SAVE PROCESSED OUTPUTS --------------------------------------------------
cat("\n=== SAVING OUTPUTS ===\n")

## 1) Seurat object (.rds)
saveRDS(so, rds_out)
cat(sprintf("Saved Seurat object   -> %s\n", rds_out))

## 2) Compact RData bundle
pca_embeddings <- Embeddings(so, "pca")  # cells x PCs
save(spatial_keep, counts_keep, hvg_genes, pca_embeddings,
     file = rdata_out)
cat(sprintf("Saved processed RData -> %s\n", rdata_out))

## 3) AnnData (.h5ad) via zellkonverter (adata.X = log-normalized)
DefaultAssay(so) <- "RNA"
mat_counts <- GetAssayData(so, assay = "RNA", layer = "counts")
mat_data   <- GetAssayData(so, assay = "RNA", layer = "data")

if (!inherits(mat_counts, "dgCMatrix")) mat_counts <- as(mat_counts, "dgCMatrix")
if (!inherits(mat_data,   "dgCMatrix")) mat_data   <- as(mat_data,   "dgCMatrix")

sce <- SingleCellExperiment(
  assays = list(
    counts    = mat_counts,
    logcounts = mat_data
  ),
  colData = S4Vectors::DataFrame(so@meta.data)
)
reducedDims(sce)$PCA <- Embeddings(so, "pca")

if (anyDuplicated(rownames(sce)) > 0) {
  sce <- sce[!duplicated(rownames(sce)), ]
}
SummarizedExperiment::rowData(sce)$highly_variable <- rownames(sce) %in% hvg_genes

if (file.exists(h5ad_out)) file.remove(h5ad_out)
zellkonverter::writeH5AD(
  sce,
  file   = h5ad_out,
  X_name = "logcounts"
)
cat(sprintf("Saved .h5ad          -> %s\n", h5ad_out))
cat("\n=== DONE ===\n")