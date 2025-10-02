# ======================== MERFISH: clean end-to-end ===========================
# Libraries
library(anndata)
library(data.table)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)

# -----------------------------------------------------------------------------
# Config
set.seed(9362)
setwd("~/Desktop/DP-RST-workload")

# -----------------------------------------------------------------------------
# Load MERFISH h5ad
adata   <- read_h5ad("./Datasets Results/Merfish_Mouse_Brain/MERFISH_0.04_20251002235120.h5ad")
X_mat   <- adata$X                            # 5488 x 155 (rows=cells, cols=genes)
obs_df  <- as.data.frame(adata$obs)           # labels: cell_class, neuron_class, domain, Region, ground_truth
gene_nm <- colnames(X_mat)
cell_id <- rownames(X_mat)

# -----------------------------------------------------------------------------
# Parse coordinates from IDs like "-3033x2825"
coord_dt <- data.table(ID = cell_id)
coord_dt[, c("x","y") := {
  s <- trimws(ID)
  i <- regexpr("x(?!.*x)", s, perl = TRUE)
  list(as.numeric(substr(s, 1, i - 1)),
       as.numeric(substr(s, i + 1, nchar(s))))
}]

# Attach x,y to obs in the same order
if (is.null(rownames(obs_df))) rownames(obs_df) <- cell_id
obs_with_xy <- obs_df %>%
  mutate(x = coord_dt$x[match(rownames(obs_df), coord_dt$ID)],
         y = coord_dt$y[match(rownames(obs_df), coord_dt$ID)])

# -----------------------------------------------------------------------------
# Build SCE (genes x cells), normalize, HVGs, PCA
X_mat <- if (!inherits(X_mat, "dgCMatrix")) as(Matrix(X_mat, sparse = TRUE), "dgCMatrix") else X_mat
counts_gc <- as(t(X_mat), "dgCMatrix")        # 155 x 5488 (rows=genes, cols=cells)

sce <- SingleCellExperiment(
  assays  = list(counts = counts_gc),
  colData = obs_with_xy,
  rowData = data.frame(gene = gene_nm, row.names = gene_nm, check.names = FALSE)
)

sce <- logNormCounts(sce)                                          # adds 'logcounts'
dec <- modelGeneVar(sce, assay.type = "logcounts")
top <- getTopHVGs(dec, n = min(2000, nrow(sce)))                   # here: all 155
sce <- runPCA(sce, subset_row = top, ncomponents = 15, exprs_values = "logcounts")

# -----------------------------------------------------------------------------
# One tidy table: PCs + coords + labels
PCA_df <- as.data.frame(reducedDim(sce, "PCA")) %>%
  mutate(ID = colnames(sce))

cd <- as.data.frame(colData(sce)) %>%
  mutate(ID = rownames(.))

pcs_labs <- PCA_df %>%
  left_join(cd[, c("ID","x","y","Region","ground_truth","cell_class","neuron_class")], by = "ID")

# -----------------------------------------------------------------------------
# Plots
p_spatial_truth <- ggplot(pcs_labs, aes(x = x, y = y, col = ground_truth)) +
  geom_point(size = 0.6) + coord_equal() +
  labs(x = "X", y = "Y", colour = "Ground truth") +
  ggtitle("MERFISH · Spatial map (ground truth)") +
  theme_minimal(base_size = 14)

p_pc_truth <- ggplot(pcs_labs, aes(x = PC1, y = PC2, col = ground_truth)) +
  geom_point(size = 0.7, alpha = 0.9) +
  labs(x = "PC1", y = "PC2", colour = "Ground truth") +
  ggtitle("MERFISH · PC1 vs PC2 (ground truth)") +
  theme_minimal(base_size = 14)

# Print
p_spatial_truth; p_pc_truth

# -----------------------------------------------------------------------------
# Save outputs
out_dir <- "./Datasets Results/Merfish_Mouse_Brain"

write.csv(pcs_labs, file.path(out_dir, "merfish_pcs_coords_labels.csv"), row.names = FALSE)
saveRDS(sce, file.path(out_dir, "merfish_sce_with_pca.rds"))

ggsave(file.path(out_dir, "spatial_ground_truth.png"), p_spatial_truth, width = 6.5, height = 5.5, dpi = 300)
# ====================== end ====================================================
