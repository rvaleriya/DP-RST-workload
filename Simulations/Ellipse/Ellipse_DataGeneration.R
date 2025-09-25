library(DP.RST)
library(sn)
library(igraph)
library(fdaPDE)
library(dplyr)
library(ggplot2)
library(plotly)
library(gridExtra)

# Set seed
set.seed(629)

setwd("~/Desktop/DP-RST-workload")

#-------------------------------------------------------------------------------
### Function to generate the boundary of the elipses
generate_ellipse <- function(center, a, b, angle = 0, n = 100) {
  t <- seq(0, 2*pi, length.out = n)
  x <- a * cos(t)
  y <- b * sin(t)

  # Apply rotation
  x_rot <- x * cos(angle) - y * sin(angle)
  y_rot <- x * sin(angle) + y * cos(angle)

  data.frame(
    x = x_rot + center[1],
    y = y_rot + center[2]
  )
}
#-------------------------------------------------------------------------------
# Parameters for the two ellipses
# Ellipse 1: centered at (0, 0), a=3, b=2, rotated by pi/6
center1 <- c(0, 0)
a1 <- 3; b1 <- 2; angle1 <- pi/6

# Ellipse 2: centered at (8, 0), a=3, b=2, rotated by -pi/8
center2 <- c(8, 0)
a2 <- 3; b2 <- 2; angle2 <- -pi/8

ellipse1 <- generate_ellipse(center1, a1, b1, angle1)
ellipse2 <- generate_ellipse(center2, a2, b2, angle2)

# Plot to verify
ggplot() +
  geom_path(data = ellipse1, aes(x = x, y = y), color = "blue") +
  geom_path(data = ellipse2, aes(x = x, y = y), color = "red") +
  coord_equal() +
  theme_minimal()

#-------------------------------------------------------------------------------
### Function to check if a point (x, y) is inside an ellipse
in_ellipse <- function(x, y, center, a, b, angle = 0) {
  # Rotate and translate coordinates to ellipse-centered system
  xp <- (x - center[1]) * cos(angle) + (y - center[2]) * sin(angle)
  yp <- -(x - center[1]) * sin(angle) + (y - center[2]) * cos(angle)
  (xp^2 / a^2 + yp^2 / b^2) <= 1
}

# Function to check if a point is at least margin 'd' away from the ellipse boundary.
safe_in_ellipse <- function(x, y, center, a, b, angle, d, n_boundary = 200) {
  boundary <- generate_ellipse(center, a, b, angle, n = n_boundary)
  distances <- sqrt((boundary$x - x)^2 + (boundary$y - y)^2)
  min(distances) >= d
}
#-------------------------------------------------------------------------------
# Generate random sample points within a bounding box that covers both ellipses.
n_points <- 1500  # increase this to get more points if needed
x_coords <- runif(n_points, min = -5, max = 12)
y_coords <- runif(n_points, min = -5, max = 5)
points <- data.frame(x = x_coords, y = y_coords)

# Filter points that are inside either ellipse
inside_points <- points[
  in_ellipse(points$x, points$y, center1, a1, b1, angle1) |
    in_ellipse(points$x, points$y, center2, a2, b2, angle2),
]

# Define the margin (minimum distance from the ellipse boundary)
margin <- 0.1

# Filter points so that they are at least 'margin' distance away from the boundary
# For each point we check the ellipse in which it lies.
points_inside <- inside_points[
  ( in_ellipse(inside_points$x, inside_points$y, center1, a1, b1, angle1) &
      mapply(safe_in_ellipse, inside_points$x, inside_points$y,
             MoreArgs = list(center = center1, a = a1, b = b1, angle = angle1, d = margin)) ) |
    ( in_ellipse(inside_points$x, inside_points$y, center2, a2, b2, angle2) &
        mapply(safe_in_ellipse, inside_points$x, inside_points$y,
               MoreArgs = list(center = center2, a = a2, b = b2, angle = angle2, d = margin)) )
  , ]


# Plot the final points/coordinates and the ellipse boundaries
ggplot() +
  geom_point(data = points_inside, aes(x = x, y = y), color = "darkgreen", alpha = 0.5) +
  geom_path(data = ellipse1, aes(x = x, y = y), color = "blue", size = 1) +
  geom_path(data = ellipse2, aes(x = x, y = y), color = "red", size = 1) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Points at Least 0.1 Units from the Ellipse Boundaries",
       x = "X Coordinate", y = "Y Coordinate")

#-------------------------------------------------------------------------------
##### GENERATE AND PLOT THE INITIAL GRAPH #####

n = nrow(points_inside)

mesh <- create.mesh.2D(points_inside)

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

graph <- trim_graph(n = nrow(points_inside), mesh = mesh, threshold = 2)

plot(graph, layout = as.matrix(points_inside), vertex.size = 1, vertex.label = NA, edge.width = 0.5)

#-------------------------------------------------------------------------------
##### OBTAIN CRUDE SPATIAL CLUSTERS #####

# Compute shortest path distances
geo_dist <- distances(graph, weights = E(graph)$weight)

### Perform clustering (e.g., hierarchical clustering)
# Convert distance matrix to "dist" object
geo_dist_dist <- as.dist(geo_dist)
geo_dist_dist[is.infinite(geo_dist_dist)] <- max(geo_dist_dist[is.finite(geo_dist_dist)], na.rm = TRUE) * 2

# Choose number of clusters
num_clusters <- 10

# Hierarchical clustering
hc <- hclust(geo_dist_dist, method = "ward.D2")

# Cut the dendrogram
clusters <- cutree(hc, k = num_clusters)

# Add clusters to coordinates
coords_geo_cluster <- points_inside %>%
  mutate(cluster = as.factor(clusters))

# Plot interactive
p <- ggplot(coords_geo_cluster, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 2) +
  theme_minimal() +
  ggtitle("Graph-Based Geodesic Clustering")

p_interactive <- ggplotly(p)
p_interactive

#-------------------------------------------------------------------------------
##### OBTAIN REFINED CLUSTERS #####

# Get unique cluster IDs
unique_clusters <- unique(coords_geo_cluster$cluster)

# Randomly assign each of the 60 clusters into 10 super-clusters (refined clusters)
shuffled <- sample(unique_clusters)
split_groups <- split(shuffled, rep(1:5, length.out = length(shuffled)))

super_cluster_map <- data.frame(
  cluster = unlist(split_groups),
  super_cluster = rep(1:5, times = sapply(split_groups, length))
)

# Merge back into the main data
coords_geo_cluster <- coords_geo_cluster %>%
  left_join(super_cluster_map, by = "cluster")

table(coords_geo_cluster$super_cluster)

# Prepare for plotting
coords_geo_cluster_plot <- coords_geo_cluster %>%
  rename(x = x, y = y)

# Plot interactive
p <- ggplot(coords_geo_cluster_plot, aes(x = x, y = y, color = as.factor(super_cluster))) +
  geom_point(size = 2) +
  theme_minimal() +
  ggtitle("Graph-Based Geodesic Clustering → 5 Super-Clusters")

p_interactive <- ggplotly(p)

p_interactive

#-------------------------------------------------------------------------------
##### GENERATE FEATURES #####

cluster_means <- matrix(c(
  0.0321, -0.399, -0.369, 0.232, -0.000643, 0.0499, -0.0663, -0.132, 0.0314, 0.0466,
  0.134, 0.595, 0.617, -0.182, -0.133, -0.301, 0.345, 0.0901, -0.158, -0.0275,
  -0.407, 0.998, 0.531, -1.49, -0.0574, 0.0806, -0.915, 0.878, 0.120, -0.599,
  1.02, 1.38, 0.375, 0.885, 0.909, 0.641, -0.798, -0.843, 0.951, 0.745,
  -1.40, 0.167, 0.654, -0.437, 0.934, 1.25, 0.325, 0.438, 0.174, 0.188
), nrow = 5, byrow = TRUE)

# Set skewness for each dimension (e.g., mostly right-skewed)
skewness_vec <- runif(10, min = 2, max = 6)  # right-skewed
skewness_vec

# [1] 3.850728 5.609784 5.975241 5.711880 3.329271 4.216154 3.483100 4.716866 4.141149
# [10] 4.453668

# Initialize feature matrix
num_spots <- nrow(coords_geo_cluster_plot)
feature_matrix <- matrix(NA, nrow = num_spots, ncol = 10)

# Generate features row-wise
for (k in 1:5) {
  idx <- which(coords_geo_cluster_plot$super_cluster == k)
  n_k <- length(idx)
  
  for (d in 1:10) {
    feature_matrix[idx, d] <- rsn(n_k, xi = cluster_means[k, d], omega = 1, alpha = skewness_vec[d])
  }
}

# Bind features to your main dataframe
colnames(feature_matrix) <- paste0("PC", 1:10)
sim_data <- cbind(coords_geo_cluster_plot, feature_matrix)
head(sim_data)

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
  sim_data,
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
ggsave("./New_Simulations/Ellipse/Ellips_DensityFeatures.pdf", 
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
  ggplot(sim_data, aes(x = x, y = y, col = .data[[var]])) +
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
ggsave("./New_Simulations/Ellipse/Ellips_SpatialFeatures.pdf", 
       plot = combined_plot, 
       width = 20, height = 5)

#-------------------------------------------------------------------------------
##### SAVE THE DATA #####

write.csv(sim_data, "./New_Simulations/Ellipse/Ellips_sim_data.csv", row.names = FALSE)

