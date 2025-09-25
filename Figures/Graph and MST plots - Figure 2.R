# Load necessary libraries
library(ggplot2)
library(sp)
library(sf)
library(gdistance)
library(rprojroot)
library(MASS)
library(ggpubr)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

### Generate random points inside a specific boundary ###

# Set seed for reproducibility
set.seed(123)  # For reproducibility

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

# Generate random points
n_points <- 200
x_random <- runif(n_points, min(x_boundary_outer), max(x_boundary_outer))
y_random <- runif(n_points, min(y_boundary_outer), max(y_boundary_outer))

# Filter points inside the boundary
inside <- is_inside(x_random, y_random, boundary)
points_inside <- data.frame(x = x_random[inside], y = y_random[inside])

# Plot the boundary and points
ggplot() +
  geom_path(data = boundary, aes(x = x, y = y), color = 'blue') +
  geom_point(data = points_inside, aes(x = x, y = y), color = 'red') +
  ggtitle("Random Points Inside the Boundary") +
  theme_minimal()


##### Points are on the grid in the Swiss Roll #####

# Set seed for reproducibility
set.seed(123) 

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
  gDistance(points_sp, boundary_sp, byid = TRUE) < 0.1 #1e-6
}

# Generate grid points within the bounding box of the boundary
x_range <- seq(floor(min(boundary$x)), ceiling(max(boundary$x)), by = 0.6)
y_range <- seq(floor(min(boundary$y)), ceiling(max(boundary$y)), by = 0.6)
grid_points <- expand.grid(x = x_range, y = y_range)

# Filter points inside the boundary
inside <- is_inside(grid_points$x, grid_points$y, boundary)
points_inside <- grid_points[inside, ]

# Remove points on the boundary
on_boundary <- is_on_boundary(points_inside$x, points_inside$y, boundary)
points_inside <- points_inside[!on_boundary, ]

# Plot the boundary and points with gradient color
points_plot <- ggplot() +
  geom_path(data = boundary, aes(x = x, y = y), color = 'navyblue') +
  geom_point(data = points_inside, aes(x = x, y = y, color = y)) +
  scale_color_gradient(low = "blue", high = "red") +
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

##### Create graph and MST #####

# Load custom functions
source("./Codes/Custom Functions/ComplexDomainFun.R")
source("./Codes/Custom Functions/plotting_fun.R")
source("./Codes/Custom Functions/BASTFun_2.R")

n = nrow(points_inside)
p = ncol(points_inside)

mesh = gen2dMesh(points_inside, boundary)
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
mstgraph = mst(graph0)  # initial spanning tree
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph, 'weight')

# plot spatial graph
spatial_graph <- plotGraph(points_inside, graph0) +
  geom_boundary(boundary, color = 'navyblue') +
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
spatial_graph

### Cluster the data points ###
kmeans_result <- kmeans(points_inside[, c("x", "y")], centers = 3)  # Change centers to the desired number of clusters
clusters <- factor(kmeans_result$cluster)

plotGraph(points_inside, graph0) +
  geom_path(data = boundary, aes(x = x, y = y), color = 'navyblue') +
  labs(x = 'Coordinate X', y = 'Coordinate Y') +
  ggtitle('Spatial Graph with Clusters') +
  geom_point(data = data.frame(points_inside, clusters), aes(x = x, y = y, color = clusters), size = 2) +
  scale_color_manual(values = c("dodgerblue", "gold1", "magenta")) +  # Manually set colors for clusters
  theme_void() +  # Remove background
  theme(legend.position = "none")


##### Run BAST DPM #####
Y_std = points_inside
# Define the "temperatures"
temp = 1
M = length(temp)

k_max = 5
mu = list() # list for initial values of mu
sigmasq_y = list()
mstgraph_lst = list()  # initial spanning trees

cluster_ids <- unique(as.numeric(clusters))
cluster = matrix(as.numeric(clusters), nrow = n, ncol = 1, byrow = F) #matrix(1, nrow = n, ncol = M)  # initial cluster memberships among 5 clusters
cluster_means_matrix <- matrix(0, nrow = length(cluster_ids), ncol = p)
# Get the means of each cluster
for (i in 1:length(cluster_ids)) {
  cluster_means_matrix[i, ] <- colMeans(Y_std[clusters == cluster_ids[i], , drop = FALSE])
}


# For each temperature, we have same initial cluster assignment
for(m in 1:M) {
  mu[[m]] = cluster_means_matrix
  sigmasq_y[[m]] = cov(Y_std) #Covariance matrix
  mstgraph_lst[[m]] = mstgraph0
}


init_val = list()
init_val[['trees']] = mstgraph_lst
init_val[['mu']] = mu 
init_val[['cluster']] = cluster #matrix(close_truth, nrow = length(close_truth), ncol = M, byrow = FALSE) #matrix(rep(kmeans(Y_std, 4)$cluster), nrow = n, ncol = M, byrow = F) #cluster #as.matrix(true_clusters) # as.matrix(kmeans(Y_std, 4)$cluster) # cluster
init_val[['sigmasq_y']] = sigmasq_y

hyperpar = list()
hyperpar[['sigmasq_mu']] = (0.5/(2*sqrt(1)))^2 #(0.5/(2*sqrt(M)))^2 # leave it to be as in univariate case (just in the code write as a matrix)
hyperpar[['lambda_s']] = diag(1, p) # lambda_s
hyperpar[['nu']] = n-p+1
hyperpar[['lambda_k']] = 8
hyperpar[['M']] = M
hyperpar[['k_max']] = k_max

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 100 # MCMC iterations
BURNIN = 0 # burnin period length
THIN = 1       # thinning intervals

dynamic_part <- "sim_test"

sim_test = fitBAST(Y_std, graph0, init_val, hyperpar, temp, 
                   MCMC, BURNIN, THIN, 
                   PT = FALSE, seed = 7394,
                   backup_d = sprintf("./Plots/%s_backup.RData", dynamic_part))


mst_only <- plotGraph(Y_std, sim_test[["tree_out"]][[1]]) +
  geom_path(data = boundary, aes(x = x, y = y), color = 'navyblue') +
  labs(x = 'Coordinate X', y = 'Coordinate Y') +
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
mst_only

mst_clust <- plotGraph(Y_std, sim_test[["tree_out"]][[1]]) +
  geom_path(data = boundary, aes(x = x, y = y), color = 'navyblue') +
  labs(x = 'Coordinate X', y = 'Coordinate Y') +
  geom_point(data = data.frame(points_inside, clusters), aes(x = x, y = y, color = clusters)) +
  scale_color_manual(values = c("dodgerblue", "gold1", "magenta")) +  # Manually set colors for clusters
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
mst_clust


combined_plots <- ggarrange(
  points_plot, spatial_graph, mst_only, mst_clust, 
  ncol = 4, nrow = 1, 
  labels = c("(a)", "(b)", "(c)", "(d)"),
  label.x = -0.03, # Adjust the horizontal position of the labels
  label.y = 1, # Adjust the vertical position of the labels
  hjust = -0.6, # Horizontal justification of the labels
  vjust = 2     # Vertical justification of the labels
)
combined_plots

ggsave("./Plots/spatial_explain_plots_paper.png", plot = combined_plots, 
       width = 40, height = 10, scale = 0.5, units = "in", dpi = 300, bg = "transparent")

