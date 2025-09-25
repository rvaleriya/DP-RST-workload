###############################################################################
# Visium-HD Mouse Embryo (square_008 µm) – end-to-end pre-processing pipeline #
###############################################################################

##### 0.  Libraries ###########################################################
library(Seurat)             # core ST object handling
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(scCustomize)        # Read10X_h5_GEO()
library(arrow)              # parquet reader
library(scran)              # HVG detection
library(scater)             # PCA
library(BayesSpace)         # spatialPreprocess()

##### 1.  Project root ########################################################
setwd("~/Desktop/DP-RST-workload/HD_Mouse_Embryo/Visium_Data")                    
message("Working dir: ", getwd())
set.seed(101)

##### 2.  Untar raw 10x archives ##############################################
# !! Run ONCE, comment out afterwards !!
# system("tar -xzf Visium_HD_Mouse_Embryo_binned_outputs.tar.gz")
# system("tar -xzf Visium_HD_Mouse_Embryo_spatial.tar.gz")

##### 3.  Define paths ########################################################
dataset_root <- getwd()
binned_dir   <- file.path(dataset_root, "binned_outputs", "square_008um")  # 8 µm bin size
spatial_dir  <- file.path(dataset_root, "binned_outputs", "square_008um", "spatial")

##### 4.  Load expression matrix ##############################################
expression_matrices <- scCustomize::Read10X_h5_GEO(data_dir = binned_dir)
counts_mat          <- expression_matrices[[1]]          # dgCMatrix (genes × spots)

##### 5.  Read spatial coordinates ###########################################
tissue_pos <- arrow::read_parquet(file.path(spatial_dir, "tissue_positions.parquet"))
tissue_pos <- as.data.frame(tissue_pos)
dim(tissue_pos) # 702244  by    6
colnames(tissue_pos) <- c("barcode", "in_tissue", "row", "col", "imagerow", "imagecol")
rownames(tissue_pos) <- tissue_pos$barcode

# ---- keep barcodes with in_tissue == 1 --------------------------------------
in_tissue_bcs <- intersect(tissue_pos$barcode[tissue_pos$in_tissue == 1],
                           colnames(counts_mat))      # safety: intersection

tissue_pos <- tissue_pos[in_tissue_bcs, ]             # reorder to match counts
dim(tissue_pos) # 344021      6
counts_mat <- counts_mat[ , in_tissue_bcs]            # subset columns

ggplot(tissue_pos, aes(x = imagecol, y = imagerow)) +
  geom_point(size = 0.15, alpha = 0.6) +
  coord_fixed() + scale_y_reverse() + theme_void() +
  ggtitle("Visium-HD Mouse Embryo – spots kept (in_tissue = 1)")

##### 6.  Basic QC – drop low-UMI spots #######################################
expr_sums     <- colSums(counts_mat)
low_cov_idx   <- which(expr_sums <= 100)
if (length(low_cov_idx) > 0) {
  counts_mat <- counts_mat[ , -low_cov_idx]
  tissue_pos <- tissue_pos[-low_cov_idx, ]            # keep metadata aligned
}
rm(expr_sums)

dim(counts_mat) # 19059 322504
dim(tissue_pos) # 322504      6

ggplot(tissue_pos, aes(x = imagecol, y = imagerow)) +
  geom_point(size = 0.15, alpha = 0.6) +
  coord_fixed() + scale_y_reverse() + theme_void() +
  ggtitle("Visium-HD Mouse Embryo – spots kept (in_tissue = 1, ≥100 UMIs)")


##### 7.  Build Seurat object #################################################
embryo_seurat <- CreateSeuratObject(counts = counts_mat,
                                    meta.data = tissue_pos)
rm(counts_mat, expression_matrices)                         # free RAM

embryo_seurat <- NormalizeData(
  embryo_seurat,
  normalization.method = "LogNormalize"  # counts per 10k, then log1p
)

##### 8.  Convert to SingleCellExperiment ####################################
sce <- as.SingleCellExperiment(embryo_seurat, assay = "RNA")
rm(embryo_seurat)


##### 10. Identify 2 000 HVGs & run PCA #######################################
dec     <- modelGeneVar(sce)
top2000 <- getTopHVGs(dec, n = 2000)
set.seed(102)
sce <- runPCA(sce, subset_row = top2000, ncomponents = 15)  # keep >3 for downstream

##### 11. Add BayesSpace metadata #############################################
sce <- spatialPreprocess(sce, platform = "Visium", skip.PCA = TRUE)

##### 12.  Extract PCA & merge with coordinates ###############################
pca_mat <- reducedDim(sce, "PCA")                          # cells × 15 PCs
dim(pca_mat) # 322504     15
pca_df  <- as.data.frame(pca_mat[, 1:15])
pca_df$ID <- rownames(pca_df)

colnames(tissue_pos)[1] <- "ID"
merged_data <- merge(pca_df, tissue_pos, by = "ID")
rownames(merged_data) <- merged_data$ID
colnames(merged_data)[c(which(colnames(merged_data) == "imagerow"),
                        which(colnames(merged_data) == "imagecol"))] <- c("x", "y")

##### 13.  Quick spatial PC heat-maps #########################################
for (pc in paste0("PC", 1:3)) {
  g <- ggplot(merged_data, aes(x = x, y = y, colour = .data[[pc]])) +
    geom_point(size = 0.4) +
    coord_fixed() + scale_y_reverse() +
    ggtitle(paste("Spatial map –", pc)) +
    theme_void()
  print(g)
}

##### 14.  Save processed objects #############################################
save(sce,          file = "/DP-RST-workload/HD_Mouse_Embryo/Processed_Data/sce_MouseEmbryo_HD.RData")
save(merged_data,  file = "/DP-RST-workload/HD_Mouse_Embryo/Processed_Data/MouseEmbryo_HD_merged_data.RData")

