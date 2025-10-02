# Load necessary libraries
library(fdaPDE)
library(ggplot2)
library(sp)
library(rgeos)
library(dplyr)
library(MASS)
library(rprojroot)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

# Get the horseshoe boundary
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes

# Connect the start and the end of the boundary
boundary_nodes <- rbind(boundary_nodes, boundary_nodes[1,])

# Convert the boundary to a data frame
boundary <- data.frame(x = boundary_nodes[, 1], y = boundary_nodes[, 2])

# Create a function to check if points are inside the boundary
is_inside <- function(x, y, boundary) {
  boundary_sp <- SpatialPolygons(list(Polygons(list(Polygon(boundary)), "boundary")))
  points_sp <- SpatialPoints(data.frame(x, y))
  !is.na(over(points_sp, boundary_sp))
}

# Create a function to check if points are on the boundary
is_on_boundary <- function(x, y, boundary) {
  boundary_sp <- SpatialLines(list(Lines(list(Line(boundary)), "boundary")))
  points_sp <- SpatialPoints(data.frame(x, y))
  gDistance(points_sp, boundary_sp, byid = TRUE) < 0.05
}

# Generate grid points within the bounding box of the boundary
x_range <- seq(floor(min(boundary$x)), ceiling(max(boundary$x)), by = 0.15)
y_range <- seq(floor(min(boundary$y)), ceiling(max(boundary$y)), by = 0.15)
grid_points <- expand.grid(x = x_range, y = y_range)

# Filter points inside the boundary
inside <- is_inside(grid_points$x, grid_points$y, boundary)
coords <- grid_points[inside, ]

# Remove points on the boundary
on_boundary <- is_on_boundary(coords$x, coords$y, boundary)
coords <- coords[!on_boundary, ]

# Plot the boundary and the points inside it
ggplot() +
  geom_polygon(data = boundary, aes(x = x, y = y), fill = NA, color = "black") +
  geom_point(data = coords, aes(x = x, y = y), color = "blue", size = 0.5) +
  coord_equal() +
  theme_minimal()

coords$cluster <- 0

plot(coords$y ~ coords$x, xlab="x", ylab="y")

### Use locator to manually assign poits to clusters ###
points <- identify(coords$x, coords$y, labels=row.names(coords))
coords[points,]$cluster <- 1

plot(coords$y ~ coords$x, col=coords$cluster, xlab="x", ylab="y")

# save(coords, boundary, file = "./Simulations/U-shape_6k/UshapeSim_coords_boundary.RData")

### Generate Y for each cluster ###

set.seed(63482)

# Define the means for each cluster's normal distribution
cluster_means <- seq(from = -2, to = 2, length.out = 6)
p = 3

# Initialize a vector to store generated Y values
Y <- matrix(0, nrow = nrow(coords), ncol = p)

# Generate Y data from a normal distribution for each cluster
for (cluster_id in unique(coords$cluster)) {
  cluster_indices <- which(coords$cluster == cluster_id)
  n_clust = length(cluster_indices)
  
  Y[cluster_indices, ] <- mvrnorm(n = n_clust, mu = rep(cluster_means[as.numeric(cluster_id)], p), 
                                  Sigma = diag(1, p))
}

sim_data <- data.frame(coords, Y)
colnames(sim_data)[(ncol(sim_data)-2):ncol(sim_data)] <- c("Y1", "Y2", "Y3")

ggplot() +
  geom_point(aes(x = x, y = y, col = cluster), data = sim_data) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") +
  ggtitle('Spatial Cluseters') +
  theme(legend.position = "none")

# save(sim_data, boundary, file = "./Simulations/U-shape_6k/UshapeSim_coords_boundary1.RData")


### Final file ###
load("./Simulations/U-shape_6k/UshapeSim_coords_boundary1.RData")

head(sim_data)
head(boundary)

ggplot() +
  geom_point(aes(x = x, y = y, col = cluster), data = sim_data) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") +
  ggtitle('Spatial Cluseters') +
  theme(legend.position = "none")
