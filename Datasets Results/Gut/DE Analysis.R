# Load necessary libraries
library(scCustomize)
library(SingleCellExperiment)
library(tidyverse)
library(stringi)
library(kableExtra)
library(ggplot2)
library(Seurat)
library(gridExtra)
library(scran)

# Set working directory to the root of the R project

set.seed(101) # Set a random seed for reproducibility

# Define the directory where the raw data is stored
data_dir <- '~/Desktop/DP-RST-workload/Gut/Space_Ranger_Data_Gut'

# Read in the expression matrices from 10X hdf5 format
expression_matrices <- Read10X_h5_GEO(data_dir = data_dir)

# Create a Seurat object from the first matrix in the list
seurat_object = CreateSeuratObject(counts = expression_matrices[[1]])

# Normalize gene expression data using log-transformation and scaling by total expression
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize")

# Retrieve all gene names from the Seurat object
all.genes <- rownames(seurat_object)

# Scale the gene expression data for all genes
seurat_object <- ScaleData(seurat_object, features = all.genes)

head(seurat_object@assays$RNA$scale.data)

# Extract the raw gene expression matrix from the Seurat object
gene_expression_matrix <- seurat_object@assays$RNA$scale.data

head(gene_expression_matrix)

dim(gene_expression_matrix) # 31053 by 3630

### Load processed data
load("~/Desktop/DP-RST-workload/Gut/swiss_roll_wt_muscle_finaltouches1.RData")

dim(swiss_roll_wt_muscle_finaltouches1) # 2568  by  7
head(swiss_roll_wt_muscle_finaltouches1)

#Extract the barcodes (rownames) of the filtered cells
subset_barcodes <- rownames(swiss_roll_wt_muscle_finaltouches1)

# Subset the gene expression matrix using the filtered barcodes
subset_gene_expression_matrix <- gene_expression_matrix[, colnames(gene_expression_matrix) %in% subset_barcodes]

# Ensure that the order of the barcodes in gene expression matches the order in subset_coords
subset_gene_expression_matrix <- subset_gene_expression_matrix[, match(subset_barcodes, colnames(subset_gene_expression_matrix))]

dim(subset_gene_expression_matrix) # 31053 by 2568
head(subset_gene_expression_matrix)

### Add clusters assignments from DP-RST output
load("~/Desktop/DP-RST-workload/Gut/DP-RST_Outputs/Gut_DP.RST_FromNewBastPT_Version2_OutputOnly.RData")

results = Gut_DP.RST_FromNewBastPT_Version2

### Choose the partition using the least squares criterion
groups_assign = results$cluster_out
teams_assign = results$teams_out
teams_number = results$j_teams_out

n = ncol(groups_assign)

# List to save co-clustering matrices
W_cum <- matrix(0, nrow = n, ncol = n)

# Get the mode of the team
tbl <- table(teams_number)
max_freq <- max(tbl)
modes <- as.numeric(names(tbl[tbl == max_freq]))

# Get the indices where teams_number equals the mode
mode_indices <- which(teams_number == modes)

# Run the loop to create co-clustering matrix for over the selected indices
for (i in 1:length(mode_indices)) {
  # Binary matrix for groups and teams
  X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
  Z <- table(sequence(length(teams_assign[[mode_indices[i]]])), teams_assign[[mode_indices[i]]])
  
  # Get the membership of observations in each team
  obs_in_teams <- X %*% Z
  
  # Compute W such that w_ij = 1 if observation i and observation j share same team
  W <- obs_in_teams %*% t(obs_in_teams)
  
  W_cum <- W_cum + W
}

mean_matrix <- W_cum/length(mode_indices)

# Calculate the vector with Frobenius norms
norm = c()

for (i in 1:length(mode_indices)) {
  # Binary matrix for groups and teams
  X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
  Z <- table(sequence(length(teams_assign[[mode_indices[i]]])), teams_assign[[mode_indices[i]]])
  
  # Get the membership of observations in each team
  obs_in_teams <- X %*% Z
  # Compute W such that w_ij = 1 if observation i and observation j share same team
  mat <- obs_in_teams %*% t(obs_in_teams)
  
  diff_matrix <- mat - mean_matrix
  norm[i] <- norm(diff_matrix, type = "F")
}

# Fins the minimum Frobenius norm index
min_index <- mode_indices[which.min(norm)]

groups_assign_out <- groups_assign[min_index, ]
teams_assign_out <- teams_assign[[min_index]]
rm(W, W_cum, mean_matrix)


# Binary membership matrix for groups and teams
X <- table(sequence(length(groups_assign_out)), groups_assign_out)
Z <- table(sequence(length(teams_assign_out)), teams_assign_out)

# Get the membership of observations in each team
obs_in_teams <- X %*% Z
# Compute W such that w_ij = 1 if observation i and observation j share same team
obs_in_teams_vec <- obs_in_teams %*% sort(unique(teams_assign_out))
length(obs_in_teams_vec) # 2568


##### Differential Expression Analysis
# Set the number of differentially expressed (DE) genes to analyze
N_DE = 5

# Convert the gene expression matrix to a numeric matrix for analysis
data = as.matrix(subset_gene_expression_matrix)

# Identify cluster markers using the scran package.
# cluster.markers = scran::findMarkers(x=data, 
#                                      pval.type = "any", 
#                                      groups=as.numeric(as.character(obs_in_teams_vec)), 
#                                      "binom")

cluster.markers = scran::findMarkers(x=data,
                                     pval.type = "some", 
                                     groups=as.numeric(as.character(obs_in_teams_vec)), 
                                     # direction="up", # Consider focusing on up-regulated markers in each cluster
                                     lfc = 1.75,     # NEW: Add a log-fold change threshold (e.g., > 25% increase)
                                     "binom")

# Convert the cluster assignment vector to numeric form
clust_VI_stable = as.numeric(as.character(obs_in_teams_vec))
# Get the rownames (gene names) of the top N_DE marker genes from the first cluster
markers      = rownames(cluster.markers[[1]][1:N_DE,])
# Combine the top marker genes' expression levels with the cluster assignments
# Transpose the matrix and bind it with the cluster information
markers_cell = cbind(t(data[markers,]), clust_VI_stable)

# Define a custom color palette for plotting, each color corresponds to a cluster
my_palette <- c("#FCDD23FF", "#F8B100FF", "#CA697CFF", "#AB74CFFF", "#48439BFF", 
  "#D9A453FF", "#E3A8A3FF", "#6A47D8FF", "#D856A7FF")

# Define a function to create a boxplot for the DE genes
Boxplot_DE <- function(markers_cell = markers_cell, N_DE){
  # Order the markers based on the cluster information
  ord_core_stable_mat = markers_cell[order(clust_VI_stable),]
  
  # Initialize the matrix for boxplot data
  BoxPlot_Data <- NULL
  
  # Loop over the number of genes (N_DE) and bind the gene expression values to BoxPlot_Data
  for (i in 1:N_DE) {
    BoxPlot_Data <- rbind(BoxPlot_Data, 
                          cbind(ord_core_stable_mat[, c(i, N_DE+1)], markers[i]))
  }
  
  # Assign column names
  colnames(BoxPlot_Data) = c("value", "cluster", "gene")
  
  # Convert to data frame and adjust types
  BoxPlot_Data <- data.frame(BoxPlot_Data)
  BoxPlot_Data$value <- as.double(BoxPlot_Data$value)
  BoxPlot_Data$cluster <- factor(BoxPlot_Data$cluster)
  
  # Plot using ggplot2
  ggplot(BoxPlot_Data, aes(y = value, x = gene, col = cluster)) + 
    scale_color_manual(values = my_palette, name = "Cluster") +
    geom_boxplot() +
    theme_bw() + 
    theme(axis.title = element_blank(), 
          axis.text.y = element_text(size = 16, face = "bold", colour = "black"), 
          axis.text.x = element_text(size = 15, angle = 45, vjust = 1, 
                                     hjust = 1, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", colour = "black"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14))
}

# Call the boxplot function to generate the plot for the top 6 DE genes across clusters
# (Left panel of Figure 5 in the main manuscript)
Plot_Boxplot_gut = Boxplot_DE(markers_cell = markers_cell, N_DE = N_DE) 
Plot_Boxplot_gut

# Save the plot
ggsave("~/Desktop/DP-RST-workload/Gut/DE_gut.png", plot = Plot_Boxplot_gut, 
       width = 8, height = 4)


### Plot the genes expressions in the intestine ###
load("~/Desktop/DP-RST-workload/Gut/swiss_roll_wt_muscle_boundary.RData")

# Genes to plot
genes_to_plot <- c("Car1", "Hmgcs2", "Cyp2c55", "Mettl7a3", "Gjb5")

# Extract expression values for all the genes and store in a list
expression_list <- lapply(genes_to_plot, function(gene) subset_gene_expression_matrix[gene, ])

# Create a combined data frame with all the gene expressions
coords_and_expr <- data.frame(
  X = swiss_roll_wt_muscle_finaltouches1$x,  # X coordinates
  Y = swiss_roll_wt_muscle_finaltouches1$y,  # Y coordinates
  Car1 = expression_list[[1]],
  Hmgcs2 = expression_list[[2]],
  Cyp2c55 = expression_list[[3]],
  Mettl7a3 = expression_list[[4]],
  Gjb5 = expression_list[[5]]
)

# Create a list to hold plots for each gene, each with its own color scale
plot_list <- lapply(genes_to_plot, function(gene) {
  ggplot(coords_and_expr, aes(x = X, y = Y, color = get(gene))) +
    geom_point(size = 1) +
    scale_color_viridis_c(option = "plasma") +  # Individual color scale for each plot
    theme_minimal() +
    theme(
      legend.title = element_blank(),  # Remove color legend title
      legend.text = element_text(size = 7),  # Make color labels smaller
      legend.key.size = unit(0.5, "cm"),  # Reduce the size of the color legend
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 20, hjust = 0.5, face = "italic")  # Centered and italic title
    ) + 
    labs(title = gene)
})

# Arrange the plots in a 2x3 grid
grid.arrange(grobs = plot_list, nrow = 2, ncol = 3)

# Create a combined plot again, but save it this time
# (Right panel of Figure 5 in the main manuscript)
g <- grid.arrange(grobs = plot_list, nrow = 3, ncol = 2)

# Save the plot
ggsave("~/Desktop/DP-RST-workload/Gut/multi_gene_expression_plot.png", 
       plot = g, width = 8, height = 10)


