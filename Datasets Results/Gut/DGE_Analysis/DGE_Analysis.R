# ===================================================================
# Load Libraries
# ===================================================================
library(Seurat)
library(scCustomize)
library(scran)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(DP.RST)
library(patchwork) # For VlnPlot annotations
library(cowplot)   # For spatial plot grid
library(dplyr)

# Set a random seed for reproducibility
set.seed(101) 

# Set the general working directory
setwd("~/Desktop/DP-RST-workload/Datasets Results/Gut")

# ===================================================================
# Load and Subset Data
# ===================================================================

# UNZIP DATA TO A TEMPORARY DIRECTORY 
zip_file <- "Space_Ranger_Data_Gut.zip"
expected_folder_name <- "Space_Ranger_Data_Gut"
temp_data_dir <- file.path(tempdir(), expected_folder_name)

cat("Unzipping", zip_file, "to temporary directory...\n")
unzip(zip_file, exdir = tempdir()) 

# Define the data_dir path to this new temporary folder
data_dir <- temp_data_dir

# Read in the expression matrices
expression_matrices <- Read10X_h5_GEO(data_dir = data_dir)

# Create a Seurat object from the first matrix
seurat_object_full <- CreateSeuratObject(counts = expression_matrices[[1]])

# Load the barcodes for the cells to keep
load("swiss_roll_wt_muscle_finaltouches1.RData")
subset_barcodes <- rownames(swiss_roll_wt_muscle_finaltouches1)

# Subset the Seurat object BEFORE normalization
seurat_subset <- subset(seurat_object_full, cells = subset_barcodes)
seurat_subset

# Clean up large object
rm(seurat_object_full, expression_matrices)
gc()

# ===================================================================
# Normalize the Subsetted Data 
# ===================================================================
# Now, normalize only the 2,568 cells
seurat_subset <- NormalizeData(seurat_subset, normalization.method = "LogNormalize")

# ===================================================================
# Load and Assign DP-RST Clusters
# ===================================================================

# Load DP-RST clustering output 
load("Results/Gut_DP.RST_FromNewBastPT_p3_Version2_OutputOnly.RData")
results = Gut_DP.RST_FromNewBastPT_Version2

# Finding the stable partition
mode_based_partition <- partition(results, method = "mode_based", batch_size = 100)

groups_assign_out <- mode_based_partition$groups_partition
teams_assign_out <- mode_based_partition$teams_partition[mode_based_partition$groups_partition]

# ===================================================================
# Add Clusters to Seurat Object
# ===================================================================

names(teams_assign_out) <- subset_barcodes

# Add this as metadata to the Seurat object
seurat_subset <- AddMetaData(
  seurat_subset,
  metadata = as.factor(teams_assign_out),
  col.name = "dp_rst_cluster"
)

# Set the "active identity" to our new clusters for DE
Idents(seurat_subset) <- "dp_rst_cluster"
cat("Clusters added to Seurat object:\n")
print(table(seurat_subset$dp_rst_cluster))

# ===================================================================
# Quality Control Metrics 
# ===================================================================

# Calculate mitochondrial percentage if not already done
seurat_subset[["percent.mt"]] <- PercentageFeatureSet(
  seurat_subset, 
  pattern = "^mt-"  # Use "^MT-" for human data
)

# Calculate standard QC metrics for ALL 5 clusters
qc_metrics <- seurat_subset@meta.data %>%
  group_by(dp_rst_cluster) %>%
  summarise(
    n_cells = n(),
    pct_total = round(100 * n() / nrow(seurat_subset@meta.data), 1),
    mean_UMIs = round(mean(nCount_RNA), 0),
    mean_genes = round(mean(nFeature_RNA), 0),
    mean_mt_pct = round(mean(percent.mt), 1)
  ) %>%
  arrange(dp_rst_cluster)

print(qc_metrics)

# Format for LaTeX
qc_metrics_latex <- qc_metrics %>%
  rename(
    Cluster = dp_rst_cluster,
    `Cells (n)` = n_cells,
    `% of Total` = pct_total,
    `Mean UMIs` = mean_UMIs,
    `Mean Genes` = mean_genes,
    `Mean MT (%)` = mean_mt_pct
  )

print(qc_metrics_latex)

# Save
write.csv(qc_metrics_latex, "DGE_Analysis/Outputs/QC_Metrics_DP-RST_5Clusters.csv", row.names = FALSE)

# ===================================================================
# Exclude Cluster 4 (due to small size) 
# ===================================================================
cat("Original object cell count:", ncol(seurat_subset), "\n")

# Create a new object that excludes all cells identified as cluster "4"
seurat_subset_filtered <- subset(
  seurat_subset,
  subset = dp_rst_cluster != "4"
)

cat("Filtered object cell count:", ncol(seurat_subset_filtered), "\n")

# Set the active identity for the filtered object
Idents(seurat_subset_filtered) <- "dp_rst_cluster"

# ===================================================================
# Differential Expression (DE) Analysis
# ===================================================================

cat("Running Seurat's FindAllMarkers on the filtered object...\n")
de_markers <- FindAllMarkers(
  seurat_subset_filtered,  
  assay = "RNA",
  slot = "data",     
  test.use = "wilcox",   
  logfc.threshold = 0.5, 
  min.pct = 0.1,        
  only.pos = TRUE       
)

# Filter for significance (FDR < 5%) BEFORE taking top genes
de_markers_sig <- de_markers %>%
  filter(p_val_adj < 0.05)

cat("Number of significant markers per cluster (FDR < 0.05):\n")
print(table(de_markers_sig$cluster))

# NOW, get the top 5 markers from the *significant* list
top5_sig_markers <- de_markers_sig %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

cat("Top 5 significant markers per cluster:\n")
print(top5_sig_markers, n = 20)

write.csv(top5_sig_markers, 
          "DGE_Analysis/Outputs/DP-RST_Top5.csv", 
          row.names = FALSE)

# ===================================================================
# FIGURE 1: Dot Plot 
# ===================================================================

# Get all top 5 genes, ordered by cluster
top5_genes <- top5_sig_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene)

fig_dotplot_horizontal <- DotPlot(
  seurat_subset_filtered,
  features = top5_genes, 
  group.by = "dp_rst_cluster",
  cols = c("lightgrey", "firebrick"), 
  dot.scale = 8,
  col.min = 0
) +
  theme_bw() +
  theme(
    # Rotate gene names on x-axis to 45 degrees
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "italic"),
    # Style cluster numbers on y-axis
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey92"),
    # Move legend to the bottom
    legend.position = "bottom" 
  ) +
  labs(
    # title = "Cluster-Specific Marker Genes",
    x = "Gene",
    y = "Cluster"
  )

print(fig_dotplot_horizontal)


ggsave("DGE_Analysis/Outputs/DotPlot_Markers_DP-RST.pdf", 
       fig_dotplot_horizontal, 
       width = 11, height = 5, bg = "white")

# ===================================================================
# FIGURE 2: Violin Plots
# ===================================================================

# Define the palette 
my_palette <- c("#FCDD23FF", "#F8B100FF", "#CA697CFF", "#AB74CFFF", "#D9A453FF")
cluster_levels <- levels(seurat_subset_filtered$dp_rst_cluster)
my_palette_filtered <- my_palette[as.numeric(cluster_levels)]


# Arrange genes by cluster
genes_ordered <- c(
  # Cluster 1 block
  "Mettl7a3", "Krt6b", "Psca", "Prss27", "Ppbp",
  # Cluster 2 block  
  "Acaa1b", "Cyp4b1", "Ugt2b5", "Defb37", "Cyp2c65",
  # Cluster 3 block
  "Mmp3", "Cd79b", "Fcer2a", "Cd19", "Cd72",
  # Cluster 5 block
  "Dkk2", "Atp1a2", "Grem1", "Gm32592", "Ift80"
)

# Create the plot (with corrected order) 
fig_violin_organized <- VlnPlot(
  seurat_subset_filtered,
  features = genes_ordered,
  pt.size = 0,
  ncol = 5,  # 5 genes per row = one row per cluster
  cols = my_palette_filtered,
  same.y.lims = FALSE  # Allow different y-scales per gene
) & # The '&' applies the theme to all subplots
  # THEME THE SUBPLOTS
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold.italic", size = 12, hjust = 0.5), # This themes the gene names
    legend.position = "bottom"
  )

print(fig_violin_organized)

ggsave("DGE_Analysis/Outputs/Violin_plot_DP-RST.pdf", fig_violin_organized, 
       width = 16, height = 12, bg = "white")

# ===================================================================
# FIGURE 3: Spatial Marker Grid 
# ===================================================================

# Define the genes to plot 
gene_pairs <- list(
  cluster1 = c("Psca", "Ppbp"),
  cluster2 = c("Acaa1b", "Cyp4b1"),
  cluster3 = c("Cd79b", "Mmp3"),
  cluster5 = c("Atp1a2", "Gm32592")
)

# Get necessary data 
lognorm_data <- GetAssayData(seurat_subset_filtered, assay = "RNA", slot = "data")

cells_to_keep <- colnames(seurat_subset_filtered)
coords_filtered <- swiss_roll_wt_muscle_finaltouches1[cells_to_keep, c("x", "y")]

# Define plotting functions 
create_plot_with_title <- function(gene) {
  expr_vals <- lognorm_data[gene, ]
  
  plot_df <- data.frame(
    X = coords_filtered$x,
    Y = coords_filtered$y,
    Expression = expr_vals
  )
  
  ggplot(plot_df, aes(x = X, y = Y, color = Expression)) +
    geom_point(size = 1.5, alpha = 0.9) +
    scale_color_viridis_c(option = "magma", name = "Expression") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 17, face = "bold.italic",
                                hjust = 0.5, margin = margin(b = 5))
    ) +
    labs(title = gene) +
    coord_fixed()
}

create_plot_clean <- function(gene) {
  expr_vals <- lognorm_data[gene, ]
  
  plot_df <- data.frame(
    X = coords_filtered$x,
    Y = coords_filtered$y,
    Expression = expr_vals
  )
  
  ggplot(plot_df, aes(x = X, y = Y, color = Expression)) +
    geom_point(size = 1.5, alpha = 0.9) +
    scale_color_viridis_c(option = "magma", name = "Expression") +
    theme_void() +
    theme(
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.5, "cm")
    ) +
    coord_fixed()
}

# Create and assemble plots 
plots <- lapply(unlist(gene_pairs), create_plot_with_title)

legend_plot <- create_plot_clean(gene_pairs$cluster1[1])
legend <- get_legend(legend_plot)

plot_grid <- plot_grid(
  plots[[1]], plots[[3]], plots[[5]], plots[[7]], # Row 1: Psca, Acaa1b, Cd79b, Atp1a2
  plots[[2]], plots[[4]], plots[[6]], plots[[8]], # Row 2: Ppbp, Cyp4b1, Mmp3, Gm32592
  ncol = 4, nrow = 2,
  labels = NULL,
  align = "hv",
  axis = "tblr"
)

cluster_labels <- ggdraw() + 
  draw_label("Cluster 1", x = 0.105, y = 0.5, size = 13, fontface = "bold", 
             color = "grey40") +
  draw_label("Cluster 2", x = 0.330, y = 0.5, size = 13, fontface = "bold", 
             color = "grey40") +
  draw_label("Cluster 3", x = 0.550, y = 0.5, size = 13, fontface = "bold", 
             color = "grey40") +
  draw_label("Cluster 5", x = 0.775, y = 0.5, size = 13, fontface = "bold", 
             color = "grey40")

fig_final <- plot_grid(
  plot_grid,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.12)
)

# Add cluster labels at the bottom
fig_final_with_labels <- plot_grid(
  fig_final,
  cluster_labels,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

print(fig_final_with_labels)

ggsave("DGE_Analysis/Outputs/Spatial_Dist_Genes_DP-RST.pdf", 
       fig_final_with_labels, 
       width = 18, height = 9, bg = "white") 

# ===================================================================
# PART 2: LOAD CLUSTERING RESULTS FROM SpaGCN and STAMarker
# ===================================================================

cat("\n=== PART 2: Loading Clustering Results ===\n")

# Load STAMarker results
stamarker <- read.csv("Results/STAMarker_Gut_joint.csv")
stamarker_clusters <- stamarker %>%
  select(barcode, stamarker_cluster = STAMarker_label)

# Load SpaGCN results
spagcn <- read.csv("~/Desktop/DP-RST-workload/Competing methods results/SpaGCN_runs/res_realdata_csv/all_clustering_results_Gut_reduced.csv")
spagcn_5clusters <- spagcn %>%
  select(barcode = X, spagcn_cluster = refined_pred_5clusters_3pcs_without_histology)

# Add other method clusters to Seurat
seurat_barcodes <- colnames(seurat_subset)

spagcn_match <- spagcn_5clusters %>% filter(barcode %in% seurat_barcodes)
spagcn_ordered <- spagcn_match[match(seurat_barcodes, spagcn_match$barcode), ]
seurat_subset$spagcn_cluster <- as.factor(spagcn_ordered$spagcn_cluster)

stamarker_match <- stamarker_clusters %>% filter(barcode %in% seurat_barcodes)
stamarker_ordered <- stamarker_match[match(seurat_barcodes, stamarker_match$barcode), ]
seurat_subset$stamarker_cluster <- as.factor(stamarker_ordered$stamarker_cluster)

# Check cluster distributions
cat("\nCluster distributions:\n")
cat("DP-RST:\n"); print(table(seurat_subset$dp_rst_cluster))
cat("SpaGCN:\n"); print(table(seurat_subset$spagcn_cluster))
cat("STAMarker:\n"); print(table(seurat_subset$stamarker_cluster))

# ===================================================================
# PART 3: FIND MARKERS FOR SpaGCN and STAMarker
# ===================================================================

cat("\n=== PART 3: Finding Markers ===\n")

# SpaGCN markers
cat("Finding SpaGCN markers...\n")
Idents(seurat_subset) <- "spagcn_cluster"
spagcn_markers <- FindAllMarkers(
  seurat_subset,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.5,
  test.use = "wilcox"
)

# STAMarker markers
cat("Finding STAMarker markers...\n")
Idents(seurat_subset) <- "stamarker_cluster"
stamarker_markers <- FindAllMarkers(
  seurat_subset,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.5,
  test.use = "wilcox"
)

# ===================================================================
# PART 4: CALCULATE MARKER QUALITY METRICS
# ===================================================================

cat("\n=== PART 4: Calculating Marker Quality ===\n")

calculate_marker_quality <- function(markers, method_name) {
  markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    summarise(
      n_significant = n(),
      mean_log2FC = mean(avg_log2FC),
      median_log2FC = median(avg_log2FC),
      max_log2FC = max(avg_log2FC),
      mean_pval = mean(-log10(p_val_adj)),
      mean_pct_diff = mean(pct.1 - pct.2)
    ) %>%
    mutate(method = method_name)
}

dprst_quality <- calculate_marker_quality(de_markers_sig, "DP-RST")
spagcn_quality <- calculate_marker_quality(spagcn_markers, "SpaGCN")
stamarker_quality <- calculate_marker_quality(stamarker_markers, "STAMarker")

all_quality <- bind_rows(dprst_quality, spagcn_quality, stamarker_quality)
print(all_quality)

method_summary <- all_quality %>%
  group_by(method) %>%
  summarise(
    total_clusters = n(),
    mean_markers_per_cluster = mean(n_significant),
    overall_mean_log2FC = mean(mean_log2FC),
    overall_mean_pval = mean(mean_pval),
    overall_mean_specificity = mean(mean_pct_diff)
  )

print(method_summary)

# ===================================================================
# PART 5: EXTRACT TOP 5 MARKERS PER CLUSTER
# ===================================================================

cat("\n=== PART 5: Extracting Top 5 Markers ===\n")

spagcn_top5 <- spagcn_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  mutate(method = "SpaGCN")

stamarker_top5 <- stamarker_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  mutate(method = "STAMarker")

cat("\nDP-RST Top 5 Markers:\n")
print(top5_sig_markers, n = 50)
cat("\nSpaGCN Top 5 Markers:\n")
print(spagcn_top5, n = 50)
cat("\nSTAMarker Top 5 Markers:\n")
print(stamarker_top5, n = 50)

write.csv(spagcn_top5, 
          "DGE_Analysis/Outputs/SpaGCN_Top5.csv", 
          row.names = FALSE)

write.csv(stamarker_top5, 
          "DGE_Analysis/Outputs/STAMarker_Top5.csv", 
          row.names = FALSE)

# ===================================================================
# PART 5.1: DOT-PLOT FOR SaGCN
# ===================================================================
# Get all top 5 genes, ordered by cluster
spagcn_top5_genes <- spagcn_top5 %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene)

fig_dotplot_spagcn <- DotPlot(
  seurat_subset,
  features = spagcn_top5_genes, 
  group.by = "spagcn_cluster",
  cols = c("lightgrey", "firebrick"), 
  dot.scale = 8,
  col.min = 0
) +
  theme_bw() +
  theme(
    # Rotate gene names on x-axis to 45 degrees
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "italic"),
    # Style cluster numbers on y-axis
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey92"),
    # Move legend to the bottom
    legend.position = "bottom" 
  ) +
  labs(
    # title = "Cluster-Specific Marker Genes",
    x = "Gene",
    y = "Cluster"
  )

print(fig_dotplot_spagcn)


ggsave("DGE_Analysis/Outputs/DotPlot_Markers_SpaGCN.pdf", fig_dotplot_spagcn, 
       width = 11, height = 5, bg = "white")

# ===================================================================
# PART 5.2: DOT-PLOT FOR STAMarker
# ===================================================================
# Get all top 5 genes, ordered by cluster
stamarker_top5_genes <- stamarker_top5 %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene)

fig_dotplot_stamarker <- DotPlot(
  seurat_subset,
  features = stamarker_top5_genes, 
  group.by = "stamarker_cluster",
  cols = c("lightgrey", "firebrick"), 
  dot.scale = 8,
  col.min = 0
) +
  theme_bw() +
  theme(
    # Rotate gene names on x-axis to 45 degrees
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "italic"),
    # Style cluster numbers on y-axis
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey92"),
    # Move legend to the bottom
    legend.position = "bottom" 
  ) +
  labs(
    # title = "Cluster-Specific Marker Genes",
    x = "Gene",
    y = "Cluster"
  )

print(fig_dotplot_stamarker)


ggsave("DGE_Analysis/Outputs/DotPlot_Markers_STAMarker.pdf", fig_dotplot_stamarker, 
       width = 11, height = 5, bg = "white")

# ===================================================================
# PART 5.1: DOT-PLOT FOR DP-RST
# ===================================================================
# Get all top 5 genes, ordered by cluster
dprst_top5_genes <- dprst_top5 %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene)

fig_dotplot_dprst <- DotPlot(
  seurat_subset,
  features = dprst_top5_genes, 
  group.by = "dp_rst_cluster",
  cols = c("lightgrey", "firebrick"), 
  dot.scale = 8,
  col.min = 0
) +
  theme_bw() +
  theme(
    # Rotate gene names on x-axis to 45 degrees
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "italic"),
    # Style cluster numbers on y-axis
    axis.text.y = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey92"),
    # Move legend to the bottom
    legend.position = "bottom" 
  ) +
  labs(
    # title = "Cluster-Specific Marker Genes",
    x = "Gene",
    y = "Cluster"
  )

print(fig_dotplot_dprst)


ggsave("DE_Analysis/Comparison_Outputs/DotPlot_Markers_DPRST.pdf", fig_dotplot_dprst, 
       width = 11, height = 5, bg = "white")


# ===================================================================
# PART 8: STATISTICAL TESTS
# ===================================================================

cat("\n=== PART 8: Statistical Tests ===\n")

kruskal_log2fc <- kruskal.test(mean_log2FC ~ method, data = all_quality)
cat("\nKruskal-Wallis test for log2FC:\n")
print(kruskal_log2fc)

dunn_log2fc <- FSA::dunnTest(mean_log2FC ~ method, data = all_quality, 
                        method = "bonferroni")
cat("\nDunn test (pairwise comparisons):\n")
print(dunn_log2fc)
