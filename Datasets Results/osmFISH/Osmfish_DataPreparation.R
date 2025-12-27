# ======================== osmFISH: clean end-to-end ===========================
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
# Load osmFISH h5ad
adata   <- read_h5ad("./Datasets Results/osmFISH/osmfish_20251122011024.h5ad")
X_mat   <- adata$X                            # 4839 x 33 (rows=cells, cols=genes)

obs_df  <- as.data.frame(adata$obs)           # labels: "ClusterName"  "ClusterID"    "Region"       "ground_truth"

spatial_coords <- adata$obsm[["spatial"]]
colnames(spatial_coords) <- c("x", "y")

gene_nm <- colnames(X_mat)
cell_id <- rownames(X_mat)

# -----------------------------------------------------------------------------
# Build SCE (genes x cells), normalize, HVGs, PCA
counts_gc <- as(t(X_mat), "dgCMatrix")       # 155 x 5488 (rows=genes, cols=cells)

sce <- SingleCellExperiment(
  assays  = list(logcounts = log1p(counts_gc)), # Direct log transformation
  colData = obs_df,
  rowData = data.frame(gene = colnames(X_mat), row.names = colnames(X_mat))
)

# Add coordinates 
colData(sce)$x <- spatial_coords[,1]
colData(sce)$y <- spatial_coords[,2]

# Variance Modeling
dec <- modelGeneVar(sce, assay.type = "logcounts")

# Note: For osmFISH (33 genes), all genes are usually HVGs. We select all of them.
top <- rownames(sce)

# Run PCA 
sce <- runPCA(sce, subset_row = top, ncomponents = 15, exprs_values = "logcounts")

# Check results
dim(reducedDim(sce, "PCA"))

# -----------------------------------------------------------------------------
# One tidy table: PCs + coords + labels
PCA_df <- as.data.frame(reducedDim(sce, "PCA")) %>%
  mutate(ID = colnames(sce))
PCA_df$x  <- sce$x   # "sce$x" is a shortcut for colData(sce)$x
PCA_df$y  <- sce$y

obs_df$ID <- rownames(obs_df)
pcs_labs  <- dplyr::left_join(PCA_df, obs_df, by = "ID")

# -----------------------------------------------------------------------------
# Plots
p_spatial_truth <- ggplot(pcs_labs, aes(x = x, y = y, col = ground_truth)) +
  geom_point(size = 0.6) + coord_equal() +
  labs(x = "X", y = "Y", colour = "Ground truth") +
  ggtitle("osmFISH · Spatial map (ground truth)") +
  theme_minimal(base_size = 14)

# Print
p_spatial_truth

# -----------------------------------------------------------------------------
# Save outputs
out_dir <- "./Datasets Results/osmFISH"

write.csv(pcs_labs, file.path(out_dir, "osmfish_pcs_coords_labels.csv"), row.names = FALSE)
saveRDS(sce, file.path(out_dir, "osmfish_sce_with_pca.rds"))

ggsave(file.path(out_dir, "osmfish_spatial_ground_truth.png"), p_spatial_truth, width = 6.5, height = 5.5, dpi = 300)

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
# [1]  -198.8543  3891.3139  6800.5180  8759.1900  9623.3101  9968.9581 10401.0181
# [8] 10717.8621 11034.7061 11092.3141 11380.3541 11870.0222 12474.9062 12964.5742
# [15] 13425.4382 13569.4582 13742.2822 13972.7142 14347.1663 14664.0103 14779.2263
# [22] 14836.8343 14980.8543 15614.5423 17602.0184 19157.4344 20338.3985 21260.1265
# [29] 21807.4025 22181.8546 22873.1506 23535.6426 24140.5266 24284.5466 24198.1346
# [36] 23910.0946 23478.0346 23218.7986 23045.9746 22930.7586 22786.7386 22873.1506
# [43] 22901.9546 22729.1306 22412.2866 22181.8546 21951.4225 21749.7945 21490.5585
# [50] 21288.9305 21116.1065 21058.4985 21058.4985 21029.6945 20943.2825 20799.2625
# [57] 20655.2425 20309.5945 20021.5545 19647.1025 19330.2584 19128.6304 19099.8264
# [64] 19157.4344 19243.8464 19531.8865 19963.9465 20338.3985 20684.0465 20885.6745
# [71] 21000.8905 21029.6945 21029.6945 21029.6945 21029.6945 21029.6945 20972.0865
# [78] 20828.0665 20540.0265 20194.3785 19819.9265 19503.0825 19330.2584 19128.6304
# [85] 19042.2184 18754.1784 18408.5304 17746.0384 16824.3104 15297.6983 13771.0862
# [92] 12359.6902 12417.2982 12014.0422 11034.7061 10141.7821  9479.2901  8932.0141
# [99]  7693.4420  6512.4780  5302.7099  4064.1379  3084.8018  1875.0338   982.1098
# [106]   377.2257  -198.8543  -256.4623  -141.2463
# 
# $y
# [1] 49846.1776 49044.4511 48128.1923 47670.0629 47440.9982 48586.3217 48929.9188
# [8] 48700.8541 47899.1276 46982.8688 46295.6746 45952.0776 45608.4805 45379.4158
# [15] 45264.8834 46066.6099 46753.8040 47440.9982 47670.0629 47097.4011 45952.0776
# [22] 44921.2864 44234.0922 43775.9628 42516.1069 41485.3157 40683.5892 40683.5892
# [29] 40569.0569 39652.7980 39080.1363 38736.5392 38278.4098 37476.6833 35987.7627
# [36] 35300.5686 35644.1656 35529.6333 35071.5038 34040.7127 33238.9862 31864.5979
# [43] 30948.3391 30490.2096 30375.6773 30833.8067 30261.1449 30032.0802 30032.0802
# [50] 30261.1449 30146.6126 29688.4832 28886.7567 27512.3684 26481.5772 25679.8507
# [57] 24992.6566 24419.9948 24419.9948 24305.4625 23961.8654 23389.2036 22358.4125
# [64] 21098.5565 20640.4271 20411.3624 20640.4271 20640.4271 20296.8301 19380.5712
# [71] 17777.1183 15715.5359 14112.0829 11706.9034  8958.1269  6552.9475  4262.3004
# [78]  3346.0415  2544.3151  2544.3151  2429.7827  2429.7827  2315.2503  1742.5886
# [85]   597.2650  -662.5909 -1120.7203 -1235.2527 -1235.2527 -1349.7850 -1349.7850
# [92] -1235.2527 17548.0535 17777.1183 18235.2477 19151.5065 19724.1683 19953.2330
# [99] 20525.8948 21098.5565 22129.3477 22816.5419 23274.6713 23961.8654 24534.5272
# [106] 25336.2537 26367.0449 27512.3684 49846.1776

# Collect and close the polygon
boundary_df <- data.frame(x = bnd_clicks$x, y = bnd_clicks$y)
boundary_closed <- rbind(boundary_df, boundary_df[1, , drop = FALSE])

# Save the boundary coordinates
bnd_csv_path <- file.path(out_dir, "osmfish_manual_boundary.csv")
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
  ggplot2::ggtitle("osmFISH · Spatial map (ground truth) with manual boundary")

p_spatial_with_bnd

# ggsave(file.path(out_dir, "spatial_ground_truth_with_boundary.png"),
#        p_spatial_with_bnd, width = 6.5, height = 5.5, dpi = 300)

# ====================== end manual boundary ===================================
