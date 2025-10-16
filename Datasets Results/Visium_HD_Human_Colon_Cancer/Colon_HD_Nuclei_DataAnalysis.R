# ==============================================================================
# Visium HD Spatial Transcriptomics Data Analysis
# ==============================================================================
# This script processes and visualizes segmented nuclei data from Visium HD
# spatial transcriptomics Human Colon Cancer dataset
# ==============================================================================

# Load required libraries
library(Matrix)
library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(sparseMatrixStats)
library(Seurat)

# Set working directory to the data folder
setwd("~/Desktop/DP-RST-workload/Datasets Results/Visium_HD_Human_Colon_Cancer")

# Function to create output directory
create_output_dir <- function() {
  output_dir <- "nuclei_data_analysis"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  return(output_dir)
}

output_dir <- create_output_dir()

# ==============================================================================
# SECTION 1: UNDERSTANDING THE DATA FILES
# ==============================================================================
cat("\n=== UNDERSTANDING YOUR VISIUM HD DATA FILES ===\n")

# File descriptions:
cat("\n1. features.tsv - Gene/Feature Information\n")
cat("   Contains information about detected genes/features including:\n")
cat("   - gene_ids: Ensembl gene IDs\n")
cat("   - feature_types: Type of feature (Gene Expression)\n")
cat("   - genome: Reference genome used (GRCh38)\n")
cat("   - n_cells_by_counts: Number of cells expressing each gene\n")
cat("   - mean_counts: Average expression across all cells\n")
cat("   - total_counts: Total counts for each gene\n")
cat("   Dimensions: 18,085 genes × 11 attributes\n")

cat("\n2. observations.tsv - Cell/Nuclei Information\n")
cat("   Contains per-cell/nuclei metrics including:\n")
cat("   - id: Unique cell identifier\n")
cat("   - n_genes_by_counts: Number of genes detected per cell\n")
cat("   - total_counts: Total UMI counts per cell\n")
cat("   - area: Nuclear area in pixels\n")
cat("   - region: Type of region (nuclei)\n")
cat("   - annotation: Cell type annotation (e.g., Neoplasm)\n")
cat("   - filtered: Quality control flag\n")
cat("   Dimensions: 305,590 nuclei × 15 attributes\n")

cat("\n3. coordinates.tsv - Spatial Coordinates\n")
cat("   Contains x,y spatial positions for each nucleus\n")
cat("   - x: X-coordinate in tissue space\n")
cat("   - y: Y-coordinate in tissue space\n")
cat("   Dimensions: 305,590 nuclei × 3 attributes\n")

cat("\n4. counts.mtx - Gene Expression Matrix\n")
cat("   Sparse matrix in Matrix Market format\n")
cat("   - Rows: Genes (18,085)\n")
cat("   - Columns: Cells/Nuclei (305,590)\n")
cat("   - Values: Expression counts (UMIs)\n")

cat("\n5. labels.tsv - Cell Type Labels\n")
cat("   Contains cell type annotations for each nucleus\n")
cat("   Dimensions: 305,590 nuclei × 2 attributes\n")

# ==============================================================================
# SECTION 2: LOADING AND PROCESSING DATA
# ==============================================================================
cat("\n\n=== LOADING DATA FILES ===\n")

zip_path <- "Nuclei_Data/segmented_nuclei.zip"
zip_contents <- unzip(zip_path, list = TRUE)
print(zip_contents)

# Load features
features <- read.table(unz(zip_path, "output/features.tsv"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Features loaded:", nrow(features), "genes\n")

# Load observations
obs <- read.table(unz(zip_path, "output/observations.tsv"),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Observations loaded:", nrow(obs), "nuclei\n")

# Load coordinates
coords <- read.table(unz(zip_path, "output/coordinates.tsv"),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Coordinates loaded:", nrow(coords), "positions\n")

# Load labels
labels <- read.table(unz(zip_path, "output/labels.tsv"),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Labels loaded:", nrow(labels), "annotations\n")

cat("\nLoading sparse count matrix from ZIP...\n")
counts <- readMM(unz(zip_path, "output/counts.mtx"))
cat("Count matrix dimensions:", nrow(counts), "genes ×", ncol(counts), "cells\n")

# Fix dimension mismatch issue
# The error suggests counts matrix has wrong dimensions
# Expected: 18085 genes × 305590 cells
# Let's check and fix if needed
if (ncol(counts) != nrow(obs)) {
  cat("\nWARNING: Dimension mismatch detected!\n")
  cat("Matrix columns:", ncol(counts), "\n")
  cat("Observation rows:", nrow(obs), "\n")
  
  # If transpose is needed:
  if (nrow(counts) == nrow(obs) && ncol(counts) == nrow(features)) {
    cat("Transposing matrix to correct orientation...\n")
    counts <- t(counts)
  }
}

# Set row and column names
rownames(counts) <- features$gene_ids
colnames(counts) <- obs$id

# Combine all metadata
metadata <- cbind(obs, coords[, c("x", "y")])
metadata$annotation <- labels$annotation

cat("\nData successfully loaded and processed!\n")
cat("Final count matrix:", nrow(counts), "genes ×", ncol(counts), "cells\n")

# ==============================================================================
# SECTION 3: DATA QUALITY ASSESSMENT
# ==============================================================================
cat("\n\n=== DATA QUALITY METRICS ===\n")

# Basic statistics
cat("\nBasic Statistics:\n")
cat("Total nuclei:", nrow(obs), "\n")
cat("Total genes:", nrow(features), "\n")
cat("Filtered nuclei:", sum(obs$filtered == "yes"), "\n")
cat("Kept nuclei:", sum(obs$filtered == "no"), "\n")
cat("Unique annotations:", paste(unique(obs$annotation), collapse = ", "), "\n")

# Gene detection statistics
gene_stats <- data.frame(
  metric = c("Median genes/cell", "Mean genes/cell", "Min genes/cell", "Max genes/cell"),
  value = c(median(obs$n_genes_by_counts), 
            mean(obs$n_genes_by_counts),
            min(obs$n_genes_by_counts),
            max(obs$n_genes_by_counts))
)
print(gene_stats)

# UMI statistics
umi_stats <- data.frame(
  metric = c("Median UMIs/cell", "Mean UMIs/cell", "Total UMIs"),
  value = c(median(obs$total_counts), 
            mean(obs$total_counts),
            sum(obs$total_counts))
)
print(umi_stats)

# ==============================================================================
# SECTION 4: QUALITY CONTROL PLOTS
# ==============================================================================
cat("\n\n=== GENERATING QUALITY CONTROL PLOTS ===\n")

# Plot 1: Distribution of genes per cell
p1 <- ggplot(obs, aes(x = n_genes_by_counts)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = median(obs$n_genes_by_counts), 
             color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = "Distribution of Detected Genes per Nucleus",
       subtitle = paste("Median:", round(median(obs$n_genes_by_counts)), "genes"),
       x = "Number of Genes",
       y = "Number of Nuclei") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(output_dir, "01_genes_per_cell_distribution.pdf"), 
       p1, width = 10, height = 6, dpi = 300)

# Plot 2: Distribution of UMI counts per cell
p2 <- ggplot(obs, aes(x = total_counts)) +
  geom_histogram(bins = 100, fill = "darkgreen", alpha = 0.7) +
  geom_vline(xintercept = median(obs$total_counts), 
             color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(labels = comma, limits = c(0, quantile(obs$total_counts, 0.99))) +
  scale_y_continuous(labels = comma) +
  labs(title = "Distribution of UMI Counts per Nucleus",
       subtitle = paste("Median:", round(median(obs$total_counts)), "UMIs"),
       x = "Total UMI Counts",
       y = "Number of Nuclei") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(output_dir, "02_umi_distribution.pdf"), 
       p2, width = 10, height = 6, dpi = 300)

# Plot 3: Genes vs UMI counts correlation
p3 <- ggplot(obs, aes(x = total_counts, y = n_genes_by_counts)) +
  geom_hex(bins = 100) +
  scale_fill_viridis(name = "Count", trans = "log10") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = "Relationship between UMI Counts and Detected Genes",
       subtitle = paste("Correlation:", round(cor(obs$total_counts, obs$n_genes_by_counts), 3)),
       x = "Total UMI Counts",
       y = "Number of Genes Detected") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(output_dir, "03_genes_vs_umi_correlation.pdf"), 
       p3, width = 10, height = 6, dpi = 300)

# Plot 4: Nuclear area distribution
p4 <- ggplot(obs, aes(x = area, fill = filtered)) +
  geom_histogram(bins = 100, alpha = 0.7) +
  scale_fill_manual(values = c("no" = "steelblue", "yes" = "red"),
                    labels = c("no" = "Kept", "yes" = "Filtered")) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = "Distribution of Nuclear Area",
       subtitle = "Colored by filtering status",
       x = "Nuclear Area (pixels)",
       y = "Number of Nuclei",
       fill = "Status") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top")

ggsave(file.path(output_dir, "04_nuclear_area_distribution.pdf"), 
       p4, width = 10, height = 6, dpi = 300)

# ==============================================================================
# SECTION 5: SPATIAL VISUALIZATION
# ==============================================================================
cat("\n=== GENERATING SPATIAL PLOTS ===\n")

# Prepare spatial data
spatial_data <- metadata[metadata$filtered == "no", ]  # Use only non-filtered cells

# Plot 5: Spatial distribution of nuclei
p5 <- ggplot(spatial_data, aes(x = x, y = y)) +
  geom_point(size = 0.1, alpha = 0.5, color = "darkblue") +
  coord_fixed() +
  labs(title = "Spatial Distribution of Nuclei",
       subtitle = paste("Total:", nrow(spatial_data), "nuclei"),
       x = "X coordinate",
       y = "Y coordinate") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(output_dir, "05_spatial_distribution.pdf"), 
       p5, width = 10, height = 10, dpi = 300)

# Plot 6: Spatial distribution colored by total counts
p6 <- ggplot(spatial_data, aes(x = x, y = y, color = log10(total_counts + 1))) +
  geom_point(size = 0.1, alpha = 0.7) +
  scale_color_viridis(name = "log10(UMI+1)") +
  coord_fixed() +
  labs(title = "Spatial Distribution of UMI Counts",
       subtitle = "Log-transformed UMI counts",
       x = "X coordinate",
       y = "Y coordinate") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")

ggsave(file.path(output_dir, "06_spatial_umi_counts.pdf"), 
       p6, width = 12, height = 10, dpi = 300)

# Plot 7: Spatial distribution colored by number of genes
p7 <- ggplot(spatial_data, aes(x = x, y = y, color = n_genes_by_counts)) +
  geom_point(size = 0.1, alpha = 0.7) +
  scale_color_viridis(name = "# Genes") +
  coord_fixed() +
  labs(title = "Spatial Distribution of Gene Detection",
       subtitle = "Number of genes detected per nucleus",
       x = "X coordinate",
       y = "Y coordinate") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right")

ggsave(file.path(output_dir, "07_spatial_gene_counts.pdf"), 
       p7, width = 12, height = 10, dpi = 300)

# Plot 8: Spatial distribution by annotation
if (length(unique(spatial_data$annotation)) > 1) {
  p8 <- ggplot(spatial_data, aes(x = x, y = y, color = annotation)) +
    geom_point(size = 0.1, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    coord_fixed() +
    labs(title = "Spatial Distribution by Cell Type Annotation",
         x = "X coordinate",
         y = "Y coordinate") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right")
  
  ggsave(file.path(output_dir, "08_spatial_annotations.pdf"), 
         p8, width = 12, height = 10, dpi = 300)
} else {
  cat("Only one annotation type found. Skipping annotation spatial plot.\n")
}

# ==============================================================================
# SECTION 6: GENE EXPRESSION ANALYSIS
# ==============================================================================
cat("\n=== ANALYZING GENE EXPRESSION PATTERNS ===\n")

# Plot 9: Gene detection rate
gene_detection <- features %>%
  mutate(detection_rate = n_cells_by_counts / nrow(obs) * 100) %>%
  arrange(desc(detection_rate))

p9 <- ggplot(gene_detection, aes(x = detection_rate)) +
  geom_histogram(bins = 100, fill = "coral", alpha = 0.7) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = "Gene Detection Rate Distribution",
       subtitle = "Percentage of cells expressing each gene",
       x = "Detection Rate (%)",
       y = "Number of Genes") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(output_dir, "09_gene_detection_rate.pdf"), 
       p9, width = 10, height = 6, dpi = 300)

# Top expressed genes
top_genes <- gene_detection %>%
  head(20) %>%
  mutate(gene = factor(gene_ids, levels = rev(gene_ids)))

p10 <- ggplot(top_genes, aes(x = detection_rate, y = gene)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  labs(title = "Top 20 Most Frequently Detected Genes",
       x = "Detection Rate (%)",
       y = "Gene") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(output_dir, "10_top_detected_genes.pdf"), 
       p10, width = 10, height = 8, dpi = 300)

# Plot 11: Mean expression vs detection rate
p11 <- ggplot(features, aes(x = n_cells_by_counts/nrow(obs)*100, 
                            y = log10(mean_counts + 1))) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Gene Expression Level vs Detection Rate",
       x = "Detection Rate (%)",
       y = "log10(Mean Expression + 1)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(output_dir, "11_expression_vs_detection.pdf"), 
       p11, width = 10, height = 6, dpi = 300)

# ==============================================================================
# SECTION 7: ADVANCED ANALYSIS - HIGHLY VARIABLE GENES
# ==============================================================================
cat("\n=== IDENTIFYING HIGHLY VARIABLE GENES ===\n")

# Calculate gene statistics for variability
# Using only non-filtered cells
kept_cells <- obs$id[obs$filtered == "no"]
counts_filtered <- counts[, as.character(kept_cells)]
counts_filtered <- as(counts_filtered, "dgCMatrix")

# Calculate mean and variance for each gene
gene_means <- rowMeans(counts_filtered)
gene_vars <- rowVars(counts_filtered)  # no dense coercion
eps <- 1e-12
gene_cv2 <- gene_vars / (pmax(gene_means, eps)^2)

# Create dataframe for highly variable genes
hvg_data <- data.frame(
  gene = rownames(counts_filtered),
  mean = gene_means,
  variance = gene_vars,
  cv2 = gene_cv2
) %>%
  filter(mean > 0.01)  # Filter very lowly expressed genes

# Fit loess to identify highly variable genes
loess_fit <- loess(log10(variance + 1) ~ log10(mean + 1), data = hvg_data, span = 0.3)
hvg_data$expected_log_variance <- predict(loess_fit, log10(hvg_data$mean + 1))
hvg_data$residual_variance <- log10(hvg_data$variance + 1) - hvg_data$expected_log_variance

# Select top highly variable genes
hvg_data <- hvg_data[order(hvg_data$residual_variance, decreasing = TRUE), ]
hvg_data$highly_variable <- FALSE
hvg_data$highly_variable[1:min(2000, nrow(hvg_data))] <- TRUE

# Plot 12: Mean-variance relationship
p12 <- ggplot(hvg_data, aes(x = log10(mean + 1), y = log10(variance + 1))) +
  geom_point(aes(color = highly_variable), alpha = 0.3, size = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
  geom_smooth(method = "loess", color = "blue", size = 1) +
  labs(title = "Mean-Variance Relationship",
       subtitle = "Identifying highly variable genes",
       x = "log10(Mean Expression + 1)",
       y = "log10(Variance + 1)",
       color = "Highly\nVariable") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(output_dir, "12_mean_variance_relationship.pdf"), 
       p12, width = 10, height = 6, dpi = 300)

# ==============================================================================
# SECTION 8: SPATIAL DISTRIBUTION OF TOP 10 HIGHLY VARIABLE GENES
# ==============================================================================
cat("\n=== PLOTTING SPATIAL DISTRIBUTION OF TOP 10 HIGHLY VARIABLE GENES ===\n")

# Get top 10 highly variable genes
top_10_hvg <- head(hvg_data[hvg_data$highly_variable, "gene"], 10)
cat("Top 10 highly variable genes:\n")
cat(paste(1:10, top_10_hvg, sep = ". "), sep = "\n")

# Prepare data for HVG spatial plots
# Create a list to store individual plots
hvg_spatial_plots <- list()

# Generate spatial plots for each of the top 10 HVGs
for (i in 1:length(top_10_hvg)) {
  gene_name <- top_10_hvg[i]
  
  # Get expression data for this gene
  gene_idx <- which(rownames(counts_filtered) == gene_name)
  
  if (length(gene_idx) > 0) {
    # Extract expression values for kept cells
    gene_expr <- as.numeric(counts_filtered[gene_idx, ])
    
    # Create dataframe for plotting
    plot_data <- spatial_data
    plot_data$expression <- gene_expr[1:nrow(spatial_data)]
    
    # Create the plot
    p_hvg <- ggplot(plot_data, aes(x = x, y = y, color = log2(expression + 1))) +
      geom_point(size = 0.05, alpha = 0.8) +
      scale_color_gradient2(low = "lightgray", mid = "yellow", high = "red", 
                            midpoint = median(log2(plot_data$expression + 1)),
                            name = "log2(Expr+1)") +
      coord_fixed() +
      labs(title = gene_name,
           subtitle = paste("HVG rank:", i)) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 8),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.key.size = unit(0.3, "cm"),
            plot.margin = margin(2, 2, 2, 2, "pt"))
    
    hvg_spatial_plots[[i]] <- p_hvg
  } else {
    cat(paste("Warning: Gene", gene_name, "not found in filtered matrix\n"))
  }
}

# Combine all HVG spatial plots into a grid
hvg_spatial_combined <- wrap_plots(hvg_spatial_plots, ncol = 3, nrow = 4) +
  plot_annotation(
    title = "Spatial Expression Patterns of Top 10 Highly Variable Genes",
    subtitle = "Each panel shows log2-transformed expression values across the tissue",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

# Save the combined HVG spatial plot
ggsave(file.path(output_dir, "13_top10_hvg_spatial_distribution.pdf"), 
       hvg_spatial_combined, width = 18, height = 20, dpi = 300)

cat("Saved spatial distribution of top 10 HVGs to: 13_top10_hvg_spatial_distribution.pdf\n")

# ==============================================================================
# SECTION 9: SUMMARY STATISTICS TABLE
# ==============================================================================
cat("\n=== CREATING SUMMARY REPORT ===\n")

nnz   <- Matrix::nnzero(counts)
total <- as.double(nrow(counts)) * as.double(ncol(counts))
sparsity <- (1 - nnz / total) * 100

# Create comprehensive summary
summary_stats <- data.frame(
  Metric = c(
    "Total number of nuclei",
    "Number of kept nuclei",
    "Number of filtered nuclei",
    "Total number of genes",
    "Median genes per nucleus",
    "Mean genes per nucleus",
    "Median UMIs per nucleus",
    "Mean UMIs per nucleus",
    "Total UMI counts",
    "Median nuclear area",
    "Tissue area covered (X range)",
    "Tissue area covered (Y range)",
    "Sparsity of expression matrix",
    "Number of highly variable genes",
    "Top highly variable gene"
  ),
  Value = c(
    format(nrow(obs), big.mark = ","),
    format(sum(obs$filtered == "no"), big.mark = ","),
    format(sum(obs$filtered == "yes"), big.mark = ","),
    format(nrow(features), big.mark = ","),
    round(median(obs$n_genes_by_counts)),
    round(mean(obs$n_genes_by_counts)),
    round(median(obs$total_counts)),
    round(mean(obs$total_counts)),
    format(sum(obs$total_counts), big.mark = ","),
    round(median(obs$area)),
    paste(round(min(coords$x)), "-", round(max(coords$x))),
    paste(round(min(coords$y)), "-", round(max(coords$y))),
    paste0(round(sparsity, 2), "%"),
    sum(hvg_data$highly_variable),
    top_10_hvg[1]
  )
)

# Save summary table
write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)

# Save HVG list
write.csv(data.frame(
  rank = 1:10,
  gene = top_10_hvg,
  mean_expression = hvg_data$mean[match(top_10_hvg, hvg_data$gene)],
  variance = hvg_data$variance[match(top_10_hvg, hvg_data$gene)],
  residual_variance = hvg_data$residual_variance[match(top_10_hvg, hvg_data$gene)]
), file.path(output_dir, "top10_hvg_statistics.csv"), row.names = FALSE)

# Print summary
cat("\n")
print(summary_stats)

# ==============================================================================
# SECTION 10: DATA EXPORT FOR DOWNSTREAM ANALYSIS
# ==============================================================================
cat("\n=== PREPARING DATA FOR DOWNSTREAM ANALYSIS ===\n")

# Create Seurat-compatible output (if needed for Seurat analysis)
# Note: Requires Seurat package

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts_filtered,
  meta.data = spatial_data,
  project = "VisiumHD"
)

# Save Seurat object
saveRDS(seurat_obj, file.path(getwd(), "Nuclei_Data/colon_hd_nuclei_seurat_object.rds"))
cat("Seurat object saved successfully!\n")

# Save processed data as R objects
save(counts_filtered, spatial_data,
     file = file.path(getwd(), "Nuclei_Data/colon_hd_nuclei_data.RData"))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nAll plots and data have been saved to:", file.path(getwd(), output_dir), "\n")
cat("\nKey outputs:\n")
cat("- 13_top10_hvg_spatial_distribution.pdf: Grid of top 10 HVG spatial patterns\n")
cat("- individual_hvg_plots/: High-resolution individual HVG plots\n")
cat("- top10_hvg_statistics.csv: Statistics for top 10 highly variable genes\n")
