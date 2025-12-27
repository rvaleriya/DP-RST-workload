library(DP.RST)
library(igraph)
library(fdaPDE)
library(sn)
library(cluster)
library(ggplot2)
library(gridExtra)
library(plotly)
library(sp)
library(sf)
library(dplyr)
library(tidyr)

# Set seed
set.seed(423)

setwd("~/Desktop/DP-RST-workload")

#-------------------------------------------------------------------------------
##### LOAD BOUNDARY #####

# Unzip the folder
zip_file <- "./Simulations/Curl/Curl_contours_csv.zip"
unzip(zip_file, exdir = "contours_csv")

# List all CSV files
csv_files <- list.files("contours_csv", pattern = "\\.csv$", full.names = TRUE)

# Read and combine all contours into one data frame
all_contours <- lapply(seq_along(csv_files), function(i) {
  df <- read.csv(csv_files[i])
  df$contour_id <- paste0("Contour_", i)  # Add contour ID
  return(df)
}) %>% bind_rows()

# Plot using ggplot2
ggplot(all_contours, aes(x = x, y = -y, group = contour_id)) +  # Negative y to match image orientation
  geom_path(size = 1) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Contours of the Curle", x = "X", y = "Y")

#-------------------------------------------------------------------------------
##### GENERATE COORDINATES WITHIN THE SHPAE #####
# Parameters
target_points <- 1000       # Desired number of points
buffer_distance <- 5       # Distance from boundary (adjust as needed)

# Read contours and convert to list of matrices
contours_list <- lapply(csv_files, function(file) {
  df <- read.csv(file)
  as.matrix(df[, c("x", "y")])
})

# Build sf object: assume largest area is outer contour
# Convert to polygons first to compute areas
poly_list <- lapply(contours_list, function(mat) st_polygon(list(mat)))

# Get areas
areas <- sapply(poly_list, st_area)
outer_index <- which.max(areas)
inner_indices <- setdiff(seq_along(poly_list), outer_index)

# Prepare list of rings: outer + holes
rings <- c(
  list(contours_list[[outer_index]]),      # Outer ring as first element
  contours_list[inner_indices]             # Holes as list of matrices
)

# Create the polygon with holes
polygon_with_holes <- st_sf(geometry = st_sfc(
  st_polygon(rings)
))

# Shrink polygon to exclude buffer area near boundary
polygon_buffered <- st_buffer(polygon_with_holes, dist = -buffer_distance)

# Dynamically compute cellsize based on area and target points
area <- st_area(polygon_buffered)
cellsize <- sqrt(as.numeric(area) / target_points)

# Generate grid points
grid <- st_make_grid(polygon_with_holes, cellsize = cellsize, what = "centers")

# Keep points inside buffered polygon
grid_in <- grid[st_within(grid, polygon_buffered, sparse = FALSE)]

# Plot to check
plot(polygon_with_holes$geometry, col = 'lightgrey', border = 'black')
plot(grid_in, add = TRUE, pch = 20, col = 'red', cex = 0.5)

# Convert to data frame for export
grid_df <- as.data.frame(st_coordinates(grid_in))
colnames(grid_df) <- c("X", "Y")
head(grid_df)
nrow(grid_df)  # Check final number of points

#-------------------------------------------------------------------------------
mesh <- create.mesh.2D(grid_df)

# Apply trimming function to cut long edges
trim_graph <- function(n, mesh, threshold = 100) {
  coords <- mesh$nodes[1:n, ]
  edge_list <- mesh$edges  # No filtering by bnd_edges, since you have no boundary edges
  
  # Compute distances
  distance <- sqrt(rowSums((coords[edge_list[, 1], ] - coords[edge_list[, 2], ])^2))
  
  # Drop long edges
  rid_drop <- distance > threshold
  edge_list <- edge_list[!rid_drop, ]
  distance <- distance[!rid_drop]
  
  # Build graph
  graph0 <- graph_from_edgelist(edge_list, directed = FALSE)
  E(graph0)$weight <- distance
  return(graph0)
}

graph <- trim_graph(n = nrow(grid_df), mesh = mesh, threshold = 30)

plot(graph, layout = as.matrix(grid_df), vertex.size = 1, vertex.label = NA, edge.width = 0.5)

#-------------------------------------------------------------------------------
##### PLOT THE INITIAL GRAPH #####

# Start with base plot
spatial_graph <- plotGraph(coords = as.matrix(grid_df), graph = graph) +
  labs(x = 'Coordinate X', y = 'Coordinate Y') +
  ggtitle('Spatial Graph')

# Add all boundaries (outer + holes)
for (boundary_mat in contours_list) {
  boundary <- list(x = boundary_mat[, 1], y = boundary_mat[, 2])
  spatial_graph <- spatial_graph + geom_boundary(boundary, colour = "black", size = 1)
}

# Print
spatial_graph

#-------------------------------------------------------------------------------
##### OBTAIN CRUDE SPATIAL CLUSTERS #####

# Compute shortest path distances
geo_dist <- distances(graph, weights = E(graph)$weight)

### Perform clustering (e.g., hierarchical clustering)
# Convert distance matrix to "dist" object
geo_dist_dist <- as.dist(geo_dist)

# Choose number of clusters
num_clusters <- 60

# Hierarchical clustering
hc <- hclust(geo_dist_dist, method = "ward.D2")

# Cut the dendrogram
clusters <- cutree(hc, k = num_clusters)

# Add clusters to coordinates
coords_geo_cluster <- grid_df %>%
  mutate(cluster = as.factor(clusters))

# Prepare boundary as data frame
boundary_df <- data.frame(x = contours_list[[outer_index]][, 1],
                          y = contours_list[[outer_index]][, 2])

coords_geo_cluster_plot <- coords_geo_cluster %>%
  rename(x = X, y = Y)

# Plot interactive
p <- ggplot(coords_geo_cluster_plot, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 2) +
  geom_path(data = boundary_df, aes(x = x, y = y), color = 'navyblue', size = 1) +
  theme_minimal() +
  ggtitle("Graph-Based Geodesic Clustering")

p_interactive <- ggplotly(p)
p_interactive

#-------------------------------------------------------------------------------
##### OBTAIN REFINED CLUSTERS #####

# Get unique cluster IDs
unique_clusters <- unique(coords_geo_cluster$cluster)

# Randomly assign each of the 60 clusters into 10 super-clusters (refined clusters)
super_cluster_map <- data.frame(
  cluster = unique_clusters,
  super_cluster = sample(1:10, length(unique_clusters), replace = TRUE)
)

# Merge back into the main data
coords_geo_cluster <- coords_geo_cluster %>%
  left_join(super_cluster_map, by = "cluster")

table(coords_geo_cluster$super_cluster)

# Prepare for plotting
coords_geo_cluster_plot <- coords_geo_cluster %>%
  rename(x = X, y = Y)

# Plot interactive
p <- ggplot(coords_geo_cluster_plot, aes(x = x, y = y, color = as.factor(super_cluster))) +
  geom_point(size = 2) +
  geom_path(data = boundary_df, aes(x = x, y = y), color = 'navyblue', size = 1) +
  theme_minimal() +
  ggtitle("Graph-Based Geodesic Clustering → 10 Super-Clusters")

p_interactive <- ggplotly(p)

p_interactive

#-------------------------------------------------------------------------------
##### GENERATE FEATURES #####

# Parameters
num_dimensions <- 10
unique_super_clusters <- unique(coords_geo_cluster$super_cluster)

# # Generate random means for each super-cluster and dimension
# super_cluster_means <- matrix(
#   runif(length(unique_super_clusters) * num_dimensions, min = -2, max = 2),
#   nrow = length(unique_super_clusters),
#   byrow = TRUE
# )
# rownames(super_cluster_means) <- unique_super_clusters
# colnames(super_cluster_means) <- paste0("PC", 1:num_dimensions)

super_cluster_means <- matrix(0, nrow = length(unique_super_clusters), ncol = num_dimensions)

# Stronger signal for PC1-PC3
super_cluster_means[, 1:3] <- matrix(
  runif(length(unique_super_clusters) * 3, min = -2, max = 2),
  nrow = length(unique_super_clusters)
)

# Weaker signal for PC4-PC10 (means closer together)
super_cluster_means[, 4:10] <- matrix(
  runif(length(unique_super_clusters) * 7, min = -1, max = 1),
  nrow = length(unique_super_clusters)
)
super_cluster_means

# [,1]       [,2]       [,3]        [,4]       [,5]       [,6]        [,7]
# [1,] -0.7698498  0.4300798  0.8965421 -0.85172903 -0.2149746  0.7701953  0.16756262
# [2,]  0.7339669 -0.3751151  0.3796656  0.45773122 -0.7396139 -0.7860417  0.09356049
# [3,]  1.4172902  0.8032535  0.9347539 -0.23251272  0.3154226 -0.9794306 -0.63867094
# [4,]  1.0630608 -0.3507444 -0.1598109  0.73070112  0.5482649 -0.9753633  0.29788819
# [5,] -1.3726295 -1.4892868 -1.3553989  0.03639748 -0.4590920  0.8893983  0.61982892
# [6,] -0.9700237 -1.2953885 -0.6524736  0.80668044 -0.1151408  0.9268289  0.38410711
# [7,]  1.6161702  0.1226605 -1.1881461  0.48894733 -0.8975577  0.6598879 -0.52500794
# [8,]  1.9688800  1.1587830  0.7605677  0.14172617  0.7775727 -0.9974281 -0.60548065
# [9,] -1.7792498  0.3490400  1.2698682 -0.24617135  0.8155762 -0.3048752 -0.62628613
# [10,]  1.5874526  0.1486751 -1.7431250 -0.99201650 -0.4021442  0.9384432 -0.12806246
# [,8]        [,9]       [,10]
# [1,] -0.5143319 -0.28044422 -0.37451595
# [2,]  0.6983443 -0.36439928 -0.53095284
# [3,] -0.4223198 -0.15617254  0.02993355
# [4,] -0.4805612 -0.72958688  0.76762113
# [5,] -0.8556260 -0.71094565 -0.10143000
# [6,]  0.3563286 -0.62034031 -0.35236001
# [7,]  0.4674007  0.03864263 -0.30300820
# [8,] -0.9343654  0.27731239 -0.38767768
# [9,]  0.2306650  0.96296599  0.80976190
# [10,] -0.7698576  0.42442185  0.45650140

# Define skewness parameter
skew_param <- runif(num_dimensions, min = -10, max = 10)  # positive skew
skew_param
# [1]  3.345068 -8.478725  3.539594 -4.109787 -8.314736 -3.754280 -8.654152  7.299265
# [9]  7.435321  9.930456

# Define random omega (scale parameter) per dimension
# omega_param <- runif(num_dimensions, min = 0.5, max = 3)
# omega_param
# # 2.4562269 0.9067784 2.1733090 1.9859214 1.6307192 1.2258828 2.3512291 0.6990043 0.6104302 2.4902774

omega_param <- c(
  runif(3, min = 1, max = 3),     # Less noise for PC1–PC3
  runif(7, min = 2, max = 5)      # More noise for PC4–PC10
)
omega_param
# [1] 2.564981 1.325423 2.338647 3.783106 3.356863 2.871059 4.221475 2.238805 2.132516
# [10] 4.388333

# Simulate data
simulated_data <- coords_geo_cluster %>%
  rowwise() %>%
  mutate(simulated_features = list({
    cluster_idx <- which(unique_super_clusters == super_cluster)
    means <- super_cluster_means[cluster_idx, ]
    sapply(1:num_dimensions, function(d) rsn(1, xi = means[d], omega = omega_param[d], alpha = skew_param[d]))
  })) %>%
  unnest_wider(simulated_features, names_sep = "_PC")

# Clean column names
colnames(simulated_data) <- gsub("simulated_features_", "", colnames(simulated_data))

head(simulated_data)

# ---------------------------------------------------------
# Overwrite PC9 and PC10 with Noise
# ---------------------------------------------------------
# Calculate number of points
n_points <- nrow(simulated_data)

# PC9: Normal distribution, Mean = 0, Variance = 1 (sd = 1)
simulated_data$PC9 <- rnorm(n_points, mean = 0, sd = 1)

# PC10: Normal distribution, Mean = 1, Variance = 2 (sd = sqrt(2))
simulated_data$PC10 <- rnorm(n_points, mean = 1, sd = sqrt(2))

head(simulated_data)

#-------------------------------------------------------------------------------
##### DENSITY PLOTS OF EACH PC PER CLUSTER #####

# Custom colors
variable_colors <- c(
  "PC1" = "blue", "PC2" = "green", "PC3" = "orange",
  "PC4" = "pink", "PC5" = "cyan", "PC6" = "yellow",
  "PC7" = "violet", "PC8" = "red", "PC9" = "purple",
  "PC10" = "darkgreen"
)

# Melt the data
df_melted <- reshape2::melt(
  simulated_data,
  id.vars = "super_cluster",
  measure.vars = paste0("PC", 1:10),
  variable.name = "Variable",
  value.name = "Value"
)

# Make sure super_cluster is a factor
df_melted$super_cluster <- as.factor(df_melted$super_cluster)

# Compute means per group for vertical lines
means_df <- df_melted %>%
  group_by(super_cluster, Variable) %>%
  summarise(mean_value = mean(Value), .groups = "drop")

# Plot
density_plots <- ggplot(df_melted, aes(x = Value, fill = Variable)) +
  geom_density(alpha = 0.6, adjust = 1.5) +
  # Add vertical mean lines
  geom_vline(data = means_df, aes(xintercept = mean_value, color = Variable),
             linetype = "dashed", size = 0.5, show.legend = FALSE) +
  facet_grid(rows = vars(super_cluster), cols = vars(Variable), scales = "free") +
  scale_fill_manual(values = variable_colors) +
  scale_color_manual(values = variable_colors) +  # Match mean lines to fill colors
  labs(
    title = "Density Plots for Each Variable by Super-Cluster",
    x = "Value",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    strip.text = element_text(size = 10)
  )

density_plots

# Save the combined plot
ggsave("./Simulations/Curl/Curl_DensityFeatures.pdf", 
       plot = density_plots, 
       width = 15, height = 10)

#-------------------------------------------------------------------------------
##### PLOT SPATIAL DISTRIBUTION OF THE FEATURES #####

# Define color scales for each variable
color_scales <- list(
  scale_color_gradient(low = "blue", high = "red"),
  scale_color_gradient(low = "green", high = "purple"),
  scale_color_gradient(low = "orange", high = "brown"),
  scale_color_gradient(low = "pink", high = "darkblue"),
  scale_color_gradient(low = "cyan", high = "magenta"),
  scale_color_gradient(low = "yellow", high = "black"),
  scale_color_gradient(low = "violet", high = "navy"),
  scale_color_gradient(low = "red", high = "gray"),
  scale_color_gradient(low = "blue", high = "gold"),
  scale_color_gradient(low = "darkgreen", high = "lightblue")
)

# Generate plots for each PC
plot_list <- lapply(seq_along(paste0("PC", 1:10)), function(i) {
  var <- paste0("PC", i)
  ggplot(simulated_data, aes(x = X, y = Y, col = .data[[var]])) +
    geom_point(size = 0.7) +  # Slightly smaller points for clarity
    scale_y_reverse() +  # Flip the Y-axis to match spatial convention
    color_scales[[i]] +  # Apply variable-specific color scale
    labs(x = 'Coordinate X', y = 'Coordinate Y', colour = var) +
    ggtitle(paste("Spatial Distribution of", var)) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
})

# Arrange all plots in a grid
combined_plot <- grid.arrange(grobs = plot_list, ncol = 5)
combined_plot

# Save the combined plot
ggsave("./Simulations/Curl/Curl_SpatialFeatures.pdf", 
       plot = combined_plot, 
       width = 20, height = 5)


#-------------------------------------------------------------------------------
##### SAVE THE DATA #####

write.csv(simulated_data, "./Simulations/Curl/Curl_sim_data.csv", row.names = FALSE)

