# ==============================================================================
# Visium HD Nuclei Data Processing Pipeline
# ------------------------------------------------------------------------------
# Steps:
#   1. Load previously saved Visium HD nuclei data
#   2. Apply quality control filters:
#        - Keep nuclei inside tissue boundary
#        - Identify and remove high-count outliers using the Kneedle method
#        - Keep genes detected in >= 3 nuclei
#   3. Seurat preprocessing:
#        - LogNormalize
#        - Select 2,000 HVGs
#        - Run PCA (15 components)
#   4. Save outputs in multiple formats:
#        - Seurat object (.rds)
#        - Processed RData bundle
#        - AnnData (.h5ad) for Python/Scanpy
#
# Outputs:
#   - Nuclei_Data/colon_hd_nuclei_seurat_object_processed.rds
#   - Nuclei_Data/colon_hd_nuclei_data_processed.RData
#   - Nuclei_Data/colon_hd_nuclei_processed.h5ad
# ==============================================================================

library(Matrix)
library(tidyverse)
library(sparseMatrixStats)
library(Seurat)
library(sf)
library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)

set.seed(42)

# ---- Paths -------------------------------------------------------------------
setwd("~/Desktop/DP-RST-workload/Datasets Results/Visium_HD_Human_Colon_Cancer")
rdata_in   <- "Nuclei_Data/colon_hd_nuclei_data.RData"                
boundary_csv <- "boundary_outer.csv"

rds_out    <- "Nuclei_Data/colon_hd_nuclei_seurat_object_processed.rds"
rdata_out  <- "Nuclei_Data/colon_hd_nuclei_data_processed.RData"

# ---- Load your saved objects -------------------------------------------------
cat("\n=== LOADING PREVIOUSLY SAVED DATA ===\n")
load(rdata_in)  

# ---- Read boundary and build polygons ---------------------------------------
cat("\n=== READING BOUNDARY & BUILDING POLYGONS ===\n")
bnd <- read.csv(boundary_csv, check.names = FALSE)

make_poly <- function(df) {
  mat <- as.matrix(df[, c("x","y")])
  if (!all(mat[1, ] == mat[nrow(mat), ])) mat <- rbind(mat, mat[1, ])
  st_polygon(list(mat))
}
poly_list <- lapply(split(bnd, bnd$polygon_id), make_poly)
polys_sfc <- st_sfc(poly_list)
polys_sf  <- st_sf(polygon_id = as.integer(names(split(bnd, bnd$polygon_id))),
                   geometry   = polys_sfc)
cat(sprintf("Constructed %d boundary polygon(s).\n", nrow(polys_sf)))

# ---- Spatial filter: inside boundary ----------------------------------------
cat("\n=== FILTERING NUCLEI INSIDE BOUNDARY ===\n")
pts_sf <- st_as_sf(spatial_data, coords = c("x","y"), remove = FALSE, crs = sf::NA_crs_)
inside <- lengths(st_intersects(pts_sf, polys_sf)) > 0
cat(sprintf("Inside boundary: %d / %d nuclei\n", sum(inside), nrow(spatial_data)))

# ---- UMI filter------------------------------------------
hist(spatial_data$total_counts, breaks = 100)
abline(v = median(spatial_data$total_counts), col = "red", lwd = 2, lty = 2)

## Kneedle (max distance to chord) for upper cutoff
tc_in <- spatial_data$total_counts[inside]
tc_sorted <- sort(tc_in, decreasing = TRUE)
n <- length(tc_sorted)
x <- seq_len(n)
xn <- (x - x[1]) / (x[n] - x[1])
yn <- (tc_sorted - tc_sorted[n]) / (tc_sorted[1] - tc_sorted[n])
num <- abs((-1) * xn - (1) * yn + 1)               # distance numerator to line (0,1)->(1,0)
den <- sqrt(2)
d <- num / den
knee_idx <- which.max(d)
upper_kneedle <- tc_sorted[knee_idx] # 813
cat(sprintf("Kneedle cutoff (upper): %.2f (idx=%d of %d)\n", upper_kneedle, knee_idx, n))

# ---- Combine filters (INSIDE & <= Kneedle) ----------------------------
keep_upper     <- spatial_data$total_counts <= upper_kneedle
keep_final     <- inside & keep_upper
cat(sprintf("Inside boundary        : %d (%.2f%%)\n", sum(inside), 100*mean(inside)))
cat(sprintf("Count <= %.2f          : %d (%.2f%%)\n", upper_kneedle, sum(keep_upper), 100*mean(keep_upper)))
cat(sprintf("Kept after ALL filters : %d (%.2f%%)\n", sum(keep_final), 100*mean(keep_final)))

# ---- Subset to kept nuclei -------------------------
kept_ids     <- spatial_data$id[keep_final]
spatial_keep <- spatial_data[keep_final, , drop = FALSE]
counts_keep  <- counts_filtered[, as.character(kept_ids), drop = FALSE]

# 1) Hard-align meta to counts order
spatial_keep <- spatial_keep[match(kept_ids, spatial_keep$id), , drop = FALSE]
stopifnot(identical(as.character(kept_ids), as.character(spatial_keep$id)))

# 2) Set barcodes as rownames for Seurat alignment
rownames(spatial_keep) <- as.character(kept_ids)

# 3) Ensure numeric x/y (THIS is the correct syntax)
spatial_keep$x <- as.numeric(spatial_keep$x)
spatial_keep$y <- as.numeric(spatial_keep$y)

# 4) Final consistency checks
stopifnot(sum(is.na(spatial_keep$x)) == 0L, sum(is.na(spatial_keep$y)) == 0L)
stopifnot(identical(colnames(counts_keep), rownames(spatial_keep)))
# ---- Gene-level QC -------------------------
cat("\n=== GENE-LEVEL QC ===\n")
cat(sprintf("Genes before filtering: %d\n", nrow(counts_keep)))
gene_detection <- rowSums(counts_keep > 0)  # number of spots where each gene is detected
min_spots <- 3  # minimum spots per gene (adjust as needed)
genes_keep <- gene_detection >= min_spots
counts_keep <- counts_keep[genes_keep, ]
cat(sprintf("Genes after filtering (detected in >=%d spots): %d\n", min_spots, nrow(counts_keep)))
cat(sprintf("Genes removed: %d (%.2f%%)\n", sum(!genes_keep), 100*mean(!genes_keep)))

# ---- Seurat preprocessing ----------------------------------------------------
cat("\n=== SEURAT PREPROCESSING (LogNormalize -> 2,000 HVGs -> 15 PCs) ===\n")
so <- CreateSeuratObject(counts = as(counts_keep, "dgCMatrix"),
                         meta.data = spatial_keep,
                         project = "VisiumHD_processed",
                         assay = "RNA",
                         min.cells = 0, min.features = 0)

so <- NormalizeData(so, normalization.method = "LogNormalize",
                    scale.factor = 1e4, verbose = FALSE)

so <- FindVariableFeatures(so, selection.method = "vst",
                           nfeatures = 2000, verbose = FALSE)
hvg_genes <- VariableFeatures(so)

so <- ScaleData(so, features = hvg_genes, verbose = FALSE)

so <- RunPCA(so, features = hvg_genes, npcs = 15, verbose = FALSE)

cat(sprintf("Computed HVGs: %d | PCs: %d\n", length(hvg_genes), 15))

# ---- SAVE PROCESSED OUTPUTS --------------------------------------------------
cat("\n=== SAVING OUTPUTS ===\n")
dir.create("Nuclei_Data", showWarnings = FALSE, recursive = TRUE)

## 1) Seurat object (.rds)
saveRDS(so, rds_out)

## 2) Compact RData bundle
pca_embeddings <- Embeddings(so, reduction = "pca")  # cells x PCs
save(spatial_keep, counts_keep, hvg_genes, pca_embeddings,
     file = rdata_out)

cat(sprintf("Saved Seurat object   -> %s\n", rds_out))
cat(sprintf("Saved processed RData -> %s\n", rdata_out))

## 3) AnnData (.h5ad) via zellkonverter
DefaultAssay(so) <- "RNA"

# Pull Seurat v5 layers
mat_counts <- GetAssayData(so, assay = "RNA", layer = "counts")
mat_data   <- GetAssayData(so, assay = "RNA", layer = "data")   # log-normalized

# Ensure sparse
if (!inherits(mat_counts, "dgCMatrix")) mat_counts <- as(mat_counts, "dgCMatrix")
if (!inherits(mat_data,   "dgCMatrix")) mat_data   <- as(mat_data,   "dgCMatrix")

# Build SCE (X will be logcounts)
sce <- SingleCellExperiment(
  assays = list(
    counts    = mat_counts,
    logcounts = mat_data
  ),
  colData = S4Vectors::DataFrame(so@meta.data)
)

# Add PCA embeddings
reducedDims(sce)$PCA <- Embeddings(so, "pca")

# Enforce unique feature names
if (anyDuplicated(rownames(sce)) > 0) {
  sce <- sce[!duplicated(rownames(sce)), ]
}

# Mark HVGs
SummarizedExperiment::rowData(sce)$highly_variable <-
  rownames(sce) %in% VariableFeatures(so)

# Write .h5ad (zellkonverter versions typically don't support overwrite=)
out_h5ad <- "Nuclei_Data/colon_hd_nuclei_processed.h5ad"
if (file.exists(out_h5ad)) file.remove(out_h5ad)

zellkonverter::writeH5AD(
  sce,
  file   = out_h5ad,
  X_name = "logcounts"  # adata.X = log-normalized matrix
)

cat(sprintf("Saved .h5ad          -> %s\n", out_h5ad))
cat("\n=== DONE ===\n")