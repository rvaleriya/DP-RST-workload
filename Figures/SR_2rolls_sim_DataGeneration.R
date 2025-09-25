# Load necessary libraries
library(ggplot2)
library(sp)
library(rgeos)
library(dplyr)
library(MASS)
library(manipulate)
library(paletteer)
library(ggpubr)
library(rprojroot)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

### Generate points inside a specific boundary ###

# Set seed for reproducibility
set.seed(6392) 

##### Create the boundary and coordinates with clusters' assignments #####
# Define the number of rolls
num_rolls <- 2  # Adjust this value to change the number of rolls

# Define the boundary for the outer spiral
theta_outer <- seq(0, num_rolls * 2 * pi, length.out = 2000)
r_outer <- 1 + 0.5 * theta_outer
x_boundary_outer <- r_outer * cos(theta_outer)
y_boundary_outer <- r_outer * sin(theta_outer)

# Define the boundary for the inner spiral (constant width)
width <- 2  # Set the constant width between spirals
r_inner <- r_outer - width
x_boundary_inner <- r_inner * cos(theta_outer)
y_boundary_inner <- r_inner * sin(theta_outer)

# Ensure the inner spiral smoothly connects with the outer spiral
connection_length <- 190
x_boundary_inner_smooth <- x_boundary_inner[-(1:connection_length)]
y_boundary_inner_smooth <- y_boundary_inner[-(1:connection_length)]

# Combine both spirals into a single boundary data frame
boundary <- data.frame(x = c(x_boundary_outer, rev(x_boundary_inner_smooth)),
                       y = c(y_boundary_outer, rev(y_boundary_inner_smooth)))

# Close the boundary by connecting the end to the start
boundary <- rbind(boundary, boundary[1, ])

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
  gDistance(points_sp, boundary_sp, byid = TRUE) < 0.15 #1e-6
}

# Generate grid points within the bounding box of the boundary
x_range <- seq(floor(min(boundary$x)), ceiling(max(boundary$x)), by = 0.7)
y_range <- seq(floor(min(boundary$y)), ceiling(max(boundary$y)), by = 0.7)
grid_points <- expand.grid(x = x_range, y = y_range)

# Filter points inside the boundary
inside <- is_inside(grid_points$x, grid_points$y, boundary)
coords <- grid_points[inside, ]

# Remove points on the boundary
on_boundary <- is_on_boundary(coords$x, coords$y, boundary)
coords <- coords[!on_boundary, ]

# Determine the roll each point belongs to based on angle
coords <- coords %>%
  mutate(angle = atan2(y, x),
         roll_segment = cut(angle, breaks = seq(-pi, pi, length.out = num_rolls + 1), labels = 1:num_rolls))

# Apply k-means clustering within each roll segment
coords <- coords %>%
  group_by(roll_segment) %>%
  mutate(cluster = kmeans(cbind(x, y), centers = 3)$cluster)


# Plot the boundary and points with clusters
points_plot <- ggplot() +
  geom_path(data = boundary, aes(x = x, y = y), color = 'navyblue') +
  geom_point(data = coords, aes(x = x, y = y, color = factor(cluster))) +
  #cale_color_manual(values = c("red", "green", "blue")) +
  labs(x = 'Coordinate X', y = 'Coordinate Y') +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(colour = "black", size = 0.5,
                             arrow = arrow(type = "closed", length = unit(0.08, "inches"))),
    plot.title = element_text(size = 36)
  )
points_plot

##### Manual processing #####
# Since there is manual work to assign coordinates to clusters, this part is not reproducible. 
# The final dataset is available to load at the end of the file.

## Reassign (manually) groups of coordinates into new clusters.
## Make spatial clusters.

## Define the means for each cluster's normal distribution ##
cluster_means <- seq(from = -10, to = 10, length.out = 3)
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

# Rename the last three columns
colnames(Y)[(ncol(Y)-2):ncol(Y)] <- c("Y1", "Y2", "Y3")
sim_data <- cbind(coords, Y)

ggplot() +
  geom_path(data = boundary, aes(x = x, y = y), color = 'navyblue') +
  geom_point(data = sim_data, aes(x = x, y = y, color = Y1)) 

# save(coords, Y, boundary, file = "./Plots/SR_2rolls_sim.RData")

# Remove all objects from the environment
rm(list = ls())
.rs.restartR()


##### Load the final data #####
load("./Plots/SR_2rolls_sim.RData")
head(coords)
