# ======================== STARmap: clean end-to-end ===========================
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
# Load STARmap h5ad
adata   <- read_h5ad("./Datasets Results/Starmap_Mouse/STARmap_20180505_BY3_1k_20251008011714.h5ad")
X_mat   <- adata$X                            # 1207 x 1020 (rows=cells, cols=genes)
obs_df  <- as.data.frame(adata$obs)           # labels: cell_class, neuron_class, domain, Region, ground_truth
gene_nm <- colnames(X_mat)
cell_id <- rownames(X_mat)

colnames(obs_df)[colnames(obs_df) %in% c("X", "Y")] <- c("x", "y")

# -----------------------------------------------------------------------------
# Build SCE (genes x cells), normalize, HVGs, PCA
X_mat <- if (!inherits(X_mat, "dgCMatrix")) as(Matrix(X_mat, sparse = TRUE), "dgCMatrix") else X_mat
counts_gc <- as(t(X_mat), "dgCMatrix")        # 155 x 5488 (rows=genes, cols=cells)

sce <- SingleCellExperiment(
  assays  = list(counts = counts_gc),
  colData = obs_df,
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

obs_df$ID <- rownames(obs_df)
pcs_labs  <- dplyr::left_join(PCA_df, obs_df, by = "ID")

# -----------------------------------------------------------------------------
# Plots
p_spatial_truth <- ggplot(pcs_labs, aes(x = x, y = y, col = ground_truth)) +
  geom_point(size = 0.6) + coord_equal() +
  labs(x = "X", y = "Y", colour = "Ground truth") +
  ggtitle("STARmap · Spatial map (ground truth)") +
  theme_minimal(base_size = 14)

p_pc_truth <- ggplot(pcs_labs, aes(x = PC1, y = PC2, col = ground_truth)) +
  geom_point(size = 0.7, alpha = 0.9) +
  labs(x = "PC1", y = "PC2", colour = "Ground truth") +
  ggtitle("STARmap · PC1 vs PC2 (ground truth)") +
  theme_minimal(base_size = 14)

# Print
p_spatial_truth; p_pc_truth

# -----------------------------------------------------------------------------
# Save outputs
out_dir <- "./Datasets Results/Starmap_Mouse"

write.csv(pcs_labs, file.path(out_dir, "starmap_pcs_coords_labels.csv"), row.names = FALSE)
saveRDS(sce, file.path(out_dir, "starmap_sce_with_pca.rds"))

ggsave(file.path(out_dir, "starmap_spatial_ground_truth.png"), p_spatial_truth, width = 6.5, height = 5.5, dpi = 300)

# ======================== Manually draw & save boundary =======================
# Show the ORIGINAL sample
plot(pcs_labs$x, pcs_labs$y,
     pch = 1, col = "grey40",
     xlab = "X", ylab = "Y")

message("Click boundary points around the sample (clockwise or counterclockwise).")
message("Press ESC (or right-click depending on OS) when finished.")

# Click points to define the boundary (locator works on the current plot)
bnd_clicks <- locator(type = "l")  # draws as you click

# $x
# [1]  -177.0557 14209.2950 14147.8149  -207.7958
# 
# $y
# [1] -6571.0104 -6571.0104   168.4225   168.4225

# Collect and close the polygon
boundary_df <- data.frame(x = bnd_clicks$x, y = bnd_clicks$y)
boundary_closed <- rbind(boundary_df, boundary_df[1, , drop = FALSE])

# Save the boundary coordinates
bnd_csv_path <- file.path(out_dir, "starmap_manual_boundary.csv")
write.csv(boundary_closed, bnd_csv_path, row.names = FALSE)
message(sprintf("Saved boundary coordinates to: %s", bnd_csv_path))

# Overlay the boundary on the original spatial plot and save a figure
p_spatial_with_bnd <- p_spatial_truth +
  ggplot2::geom_path(
    data = boundary_closed,
    aes(x = x, y = y),
    linewidth = 0.9, color = "black"
  ) +
  ggplot2::geom_point(
    data = boundary_closed,
    aes(x = x, y = y),
    size = 1.2, color = "black"
  ) +
  ggplot2::ggtitle("STARmap · Spatial map (ground truth) with manual boundary")

p_spatial_with_bnd

# ggsave(file.path(out_dir, "spatial_ground_truth_with_boundary.png"),
#        p_spatial_with_bnd, width = 6.5, height = 5.5, dpi = 300)

# ====================== end manual boundary ===================================
