# Load necessary libraries
library(ggplot2)
library(sp)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(plotly)

#-------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/DP-RST-workload")
getwd()

# Load datasets
gut_df_wt_muscle <- readRDS("./Gut/gut_df_wt_muscle.rds")
load("./Gut/swiss_roll_wt_muscle_boundary.RData")

#-------------------------------------------------------------------------------

# Extract first three principal components and spatial coordinates
Y_sample <- gut_df_wt_muscle[, 1:10]
loc <- gut_df_wt_muscle[, c("x", "y")]
bnd <- boundary

# Standardize coordinates and Y values
coords <- as.data.frame(apply(loc, 2, scale))
Y_std <- as.data.frame(apply(Y_sample, 2, scale))

# Standardize the boundary given the params of coordinates
bnd_scaled <- data.frame(
  x = (bnd$x - mean(loc[, 1])) / sd(loc[, 1]),
  y = (bnd$y - mean(loc[, 2])) / sd(loc[, 2])
)

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std)

#-------------------------------------------------------------------------------
# Get the palette from paletteer
palette_colors <- paletteer::paletteer_c("grDevices::Blue-Red 2", n = 100)

# Convert the palette to a format usable by Plotly
palette_colors <- as.character(palette_colors)

#-------------------------------------------------------------------------------

# Create a colorscale that Plotly can use
colorscale <- list()
for (i in seq_along(palette_colors)) {
  colorscale[[i]] <- list((i - 1) / (length(palette_colors) - 1), palette_colors[i])
}

# Create the plot using plotly
p <- plot_ly() %>%
  add_trace(
    data = df_subset,
    x = ~(-x), z = ~y, y = 0,
    type = 'scatter3d', mode = 'markers',
    marker = list(size = 3, color = ~PC3, colorscale = colorscale, opacity = 1, 
                  colorbar = list(len = 0.8, title = "PCA values")),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = df_subset,
    x = ~(-x), z = ~y, y = 5,
    type = 'scatter3d', mode = 'markers',
    marker = list(size = 3, color = ~PC2, colorscale = colorscale, opacity = 1),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = df_subset,
    x = ~(-x), z = ~y, y = 8,
    type = 'scatter3d', mode = 'markers',
    marker = list(size = 3, color = ~PC1, colorscale = colorscale, opacity = 1),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = bnd_scaled,
    x = ~(-x), z = ~y, y = rep(0, nrow(bnd_scaled)),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'black', width = 2),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = bnd_scaled,
    x = ~(-x), z = ~y, y = rep(5, nrow(bnd_scaled)),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'black', width = 2),
    showlegend = FALSE  # Hide legend
  ) %>%
  add_trace(
    data = bnd_scaled,
    x = ~(-x), z = ~y, y = rep(8, nrow(bnd_scaled)),
    type = 'scatter3d', mode = 'lines',
    line = list(color = 'black', width = 2),
    showlegend = FALSE  # Hide legend
  ) %>%
  layout(scene = list(
    xaxis = list(title = '', showticklabels = FALSE, showbackground = FALSE, zeroline = FALSE, showgrid = FALSE),
    yaxis = list(title = '', showticklabels = FALSE, showbackground = FALSE, zeroline = FALSE, showgrid = FALSE),
    zaxis = list(title = '', showticklabels = FALSE, showbackground = FALSE, zeroline = FALSE, showgrid = FALSE),
    bgcolor = 'rgba(0,0,0,0)'
  ),
  paper_bgcolor = 'rgba(0,0,0,0)',
  plot_bgcolor = 'rgba(0,0,0,0)'
  )

p

#-------------------------------------------------------------------------------

# Compute quantile-based limits for better contrast
all_values <- unlist(Y_std)
color_limits <- range(unlist(Y_std), na.rm = TRUE)

# Function to create individual plot
create_plot <- function(pc_column, pc_name) {
  ggplot(df_subset, aes(x = x, y = y, color = .data[[pc_column]])) +
    geom_point(size = 0.5) +
    geom_path(data = bnd_scaled, aes(x = x, y = y), color = "black", linewidth = 0.3, inherit.aes = FALSE) +
    scale_color_gradientn(colors = palette_colors, limits = color_limits, oob = scales::squish) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    ggtitle(pc_name)
}

# Create plots for PC1 to PC10
plot_list <- lapply(1:3, function(i) {
  create_plot(paste0("PC", i), paste0("PC", i))
})

# Arrange plots in 2x5 grid, shared legend
final_plot <- wrap_plots(plot_list, nrow = 1, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

# Print the final plot
print(final_plot)

# Save the plot
ggsave("./Figures/PCs_figures/PCs_expression_plot_3PC_CommonColorScale.png", plot = final_plot, width = 8, height = 3, dpi = 300)


#-------------------------------------------------------------------------------

# Compute quantile-based limits for better contrast
all_values <- unlist(Y_std[, 1:3])
lower_limit <- quantile(all_values, 0.01, na.rm = TRUE)
upper_limit <- quantile(all_values, 0.99, na.rm = TRUE)
color_limits <- c(lower_limit, upper_limit)

# Function to create individual plot
create_plot <- function(pc_column, pc_name) {
  ggplot(df_subset, aes(x = x, y = y, color = .data[[pc_column]])) +
    geom_point(size = 0.5) +
    geom_path(data = bnd_scaled, aes(x = x, y = y), color = "black", linewidth = 0.3, inherit.aes = FALSE) +
    scale_color_gradientn(colors = palette_colors, limits = color_limits, oob = scales::squish) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    ggtitle(pc_name)
}

# Create plots for PC1 to PC10
plot_list <- lapply(1:10, function(i) {
  create_plot(paste0("PC", i), paste0("PC", i))
})

# Arrange plots in 2x5 grid, shared legend
final_plot <- wrap_plots(plot_list, nrow = 2, ncol = 5, guides = "collect") &
  theme(legend.position = "bottom")

# Print the final plot
print(final_plot)

# Save the plot
ggsave("./Figures/PCs_figures/PCs_expression_plot.png", plot = final_plot, width = 12, height = 6, dpi = 300)
# ggsave("./Figures/PCs_figures/PCs_expression_plot.pdf", plot = final_plot, width = 12, height = 6)

#-------------------------------------------------------------------------------
##### DENSITY PLOTS OF PCs FOR TRUE CLUSTERS #####

# Prepare data for plotting
plot_data <- gut_df_wt_muscle %>%
  select(z, PC1:PC10) %>%
  pivot_longer(cols = PC1:PC10, names_to = "PC", values_to = "Value")

# Properly rename factor levels to match desired labels
plot_data$z <- forcats::fct_recode(plot_data$z,
                                   "Edema" = "Edema phenotype",
                                   "Inflammation" = "Inflammation phenotype",
                                   "Hyperchromasy" = "Loss of goblet cells",
                                   "Lymphoid patch" = "Lymphoid phenotype",
                                   "Normal tissue" = "Normal phenotype"
)

# Reorder factor levels for cluster rows
plot_data$z <- factor(plot_data$z, levels = c(
  "Normal tissue",
  "Inflammation",
  "Edema",
  "Hyperchromasy",
  "Lymphoid patch"
))

# Soft fill colors
fill_colors <- c(
  "Normal tissue" = "#bbc6dd",
  "Inflammation" = "#f4e6aa",
  "Edema" = "#b0d99c",
  "Hyperchromasy" = "#eba8a7",
  "Lymphoid patch" = "#bab4d8"
)

# Darker outline colors
line_colors <- c(
  "Normal tissue" = "#3C4F8A",
  "Inflammation" = "#ECD900",
  "Edema" = "#5A913F",
  "Hyperchromasy" = "#C3423F",
  "Lymphoid patch" = "#8B5FBF"
)

# Order PCs correctly
plot_data$PC <- factor(plot_data$PC, levels = paste0("PC", 1:10))

# Create the density plot
PCs_densities_truth <- ggplot(plot_data, aes(x = Value, fill = z, color = z)) +
  geom_density(alpha = 0.6, adjust = 1.2, linewidth = 0.7) +
  facet_grid(rows = vars(z), cols = vars(PC), scales = "free", switch = "y") +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = line_colors) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12, face = "bold"),
    strip.text.x = element_text(size = 10),
    strip.background = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines"),
    axis.title.y = element_blank(),
    axis.title.y.right = element_blank()  # <- remove any y-axis label completely
  ) +
  labs(x = "PC value", y = NULL)

PCs_densities_truth

ggsave("./Figures/PCs_figures/PCs_densities_truth.png", plot = PCs_densities_truth, width = 12, height = 6, dpi = 300)

#-------------------------------------------------------------------------------

# Prepare data with standardized PC scores
plot_data_std <- data.frame(z = gut_df_wt_muscle$z, Y_std)

# Clean cluster names
plot_data_std$z <- forcats::fct_recode(plot_data_std$z,
                                       "Edema" = "Edema phenotype",
                                       "Inflammation" = "Inflammation phenotype",
                                       "Hyperchromasy" = "Loss of goblet cells",
                                       "Lymphoid patch" = "Lymphoid phenotype",
                                       "Normal tissue" = "Normal phenotype"
)

# Reorder cluster factor levels
plot_data_std$z <- factor(plot_data_std$z, levels = c(
  "Normal tissue",
  "Inflammation",
  "Edema",
  "Hyperchromasy",
  "Lymphoid patch"
))

# Compute mean of standardized PC scores by cluster
mean_table_std <- plot_data_std %>%
  group_by(z) %>%
  summarise(across(starts_with("PC"), mean, na.rm = TRUE))

# Print the table
print(mean_table_std)

# z                  PC1    PC2    PC3    PC4       PC5     PC6     PC7     PC8     PC9    PC10
# 
# 1 Normal tissue   0.0321 -0.399 -0.369  0.232 -0.000643  0.0499 -0.0663 -0.132   0.0314  0.0466
# 2 Inflammation    0.134   0.595  0.617 -0.182 -0.133    -0.301   0.345   0.0901 -0.158  -0.0275
# 3 Edema          -0.407   0.998  0.531 -1.49  -0.0574    0.0806 -0.915   0.878   0.120  -0.599 
# 4 Hyperchromasy   1.02    1.38   0.375  0.885  0.909     0.641  -0.798  -0.843   0.951   0.745 
# 5 Lymphoid patch -1.40    0.167  0.654 -0.437  0.934     1.25    0.325   0.438   0.174   0.188 

#-------------------------------------------------------------------------------
##### SPATIAL MAPS OF TOP GENES FROM THE PCA LOADINGS #####

# Load PCA loadings
loadings_df <- read.csv("./Gut/PCA_loadings.csv", row.names = 1)
loadings_matrix <- as.matrix(loadings_df)

# Extract top 5 genes per PC
top_genes_list <- lapply(1:3, function(pc_num) {
  pc_name <- paste0("PC", pc_num)
  top_genes <- names(sort(abs(loadings_matrix[, pc_name]), decreasing = TRUE))[1:5]
  return(top_genes)
})
names(top_genes_list) <- paste0("PC", 1:3)

# Define PCs
pcs <- names(top_genes_list)


library(scCustomize)          
library(SingleCellExperiment) 

# Prepare gene expression data
selected_genes <- unique(unlist(top_genes_list))
data_dir <- './Gut/Space_Ranger_Data_Gut'
expression_matrices <- Read10X_h5_GEO(data_dir = data_dir)
seurat_object <- CreateSeuratObject(counts = expression_matrices[[1]])
sce <- as.SingleCellExperiment(seurat_object)

# Log-transform
log_count <- log(sce@assays@data@listData$counts + 1)
M1 <- as(log_count, "CsparseMatrix")
sce@assays@data@listData$logcounts <- M1
rm(log_count)

# Extract expression
gene_expression <- as.matrix(sce@assays@data$logcounts[selected_genes, ])
gene_expression_df <- as.data.frame(t(gene_expression))
gene_expression_df$ID <- rownames(gene_expression_df)

# Merge with coordinates
df_subset$ID <- rownames(gut_df_wt_muscle)
plot_data_genes <- left_join(df_subset, gene_expression_df, by = "ID")

# Define PCs
pcs <- names(top_genes_list)

# Common color scale across all plots
color_limits <- range(unlist(gene_expression_df[selected_genes]), na.rm = TRUE)

# Generate all plots in matrix form: list of lists (rows = PCs, columns = genes)
plot_matrix <- lapply(names(top_genes_list), function(pc) {
  lapply(top_genes_list[[pc]], function(gene) {
    ggplot(plot_data_genes, aes(x = x, y = y, color = .data[[gene]])) +
      geom_point(size = 0.5) +
      geom_path(data = bnd_scaled, aes(x = x, y = y), color = "black", linewidth = 0.3, inherit.aes = FALSE) +
      scale_color_gradientn(colors = palette_colors, limits = color_limits, oob = scales::squish) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 9),
        plot.margin = margin(0, 0, 0, 0)
      ) +
      ggtitle(gene)
  })
})

# Compute global color scale limits
color_limits <- range(plot_data_genes %>% select(all_of(selected_genes)), na.rm = TRUE)

# Prepare plot matrix
plot_matrix <- lapply(pcs, function(pc) {
  plots <- list()
  
  # Add header for the row
  header_plot <- ggplot() + 
    theme_void() + 
    annotate("text", x = 0.5, y = 0.5, label = pc, size = 5, fontface = "bold") +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  plots[[1]] <- header_plot
  
  # Add gene plots
  for (gene in top_genes_list[[pc]]) {
    p <- ggplot(plot_data_genes, aes(x = x, y = y, color = .data[[gene]])) +
      geom_point(size = 0.5) +
      geom_path(data = bnd_scaled, aes(x = x, y = y), color = "black", linewidth = 0.3, inherit.aes = FALSE) +
      scale_color_gradientn(colors = palette_colors, limits = color_limits, oob = scales::squish) +
      theme_void() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        plot.margin = margin(0, 0, 0, 0)
      ) +
      ggtitle(gene)
    
    plots[[length(plots) + 1]] <- p
  }
  
  return(plots)
})


plot_list <- unlist(plot_matrix, recursive = FALSE)

final_plot <- wrap_plots(plot_list, ncol = length(top_genes_list[[1]]) + 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(final_plot)

ggsave("./Figures/PCs_figures/spatial_top5_genes_by_PC_rows.png", final_plot, width = 12, height = 6, dpi = 300)

