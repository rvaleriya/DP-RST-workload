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

# Load custom functions
source("./Codes/Custom Functions/ComplexDomainFun.R")
source("./Codes/Custom Functions/plotting_fun.R")
source("./Codes/Custom Functions/BASTFun_2.R")

# Load boundary and coordinates data
load("./Plots/SR_2rolls_sim.RData")

sim_data <- cbind(coords, Y)

# Extract first three PCAs
Y_sample = sim_data[, 4:6]
loc = sim_data[, 1:2]
bnd = boundary

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

# Standardize boundary coordinates for consistent plotting
bnd_scaled = list()
bnd_scaled$x <- (bnd$x - mean(loc[,1]))/sd(loc[,1])
bnd_scaled$y <- (bnd$y - mean(loc[,2]))/sd(loc[,2])

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std)

# Plot the boundary and points colored by the first PCA (Y1)
ggplot() +
  geom_boundary(bnd_scaled) +
  geom_point(aes(x = x, y = y, col = Y1), data = df_subset) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "Y1") +
  ggtitle('First PCA')


# Generate a 2D mesh for triangulation
mesh = gen2dMesh(coords, bnd_scaled)
graph0 = constrainedDentri(n, mesh) # Generate a constrained Delaunay triangulation graph
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id

# Create the Minimum Spanning Tree (MST) from the graph
mstgraph = mst(graph0)  # Generate the MST from the graph
graph0 = delete_edge_attr(graph0, 'weight') # Remove edge weights from graph
mstgraph0 = delete_edge_attr(mstgraph, 'weight') # Remove edge weights from the MST

# plot spatial graph
plotGraph(coords, graph0) +
  geom_boundary(bnd_scaled) +
  labs(x = 'Coordinate X', y = 'Coordinate Y') +
  ggtitle('Spatial Graph')

# Create MST for the current cluster assignments
subgraphs = list()
eid_btw_mst = list()

cluster_m = sim_data$spatial_cluster # Cluster labels from the simulation data
mstgraph_m = mstgraph0 # Initialize the MST graph for cluster assignments
clust_vid_m = split(1:n, cluster_m) # Split vertex IDs by clusters
inc_mat = get.edgelist(graph0, names = F)  # Incidence matrix (edges) from the graph
adj_list = lapply(as_adj_list(graph0), FUN = function(x) {x$vid}) # Adjacency list for vertices
adj_edge_list = lapply(as_adj_edge_list(graph0), FUN = function(x) {x$eid}) # Adjacency list for edges

# Create subgraphs for each cluster
subgraphs = lapply(clust_vid_m, function(vids, mstgraph) {
  induced_subgraph(mstgraph, vids)
}, mstgraph_m)

# Determine edge status and propose a new MST
edge_status = apply(as.matrix(cluster_m), 2, FUN = getEdgeStatus, inc_mat)
mstgraph_m = proposeMST(graph0, edge_status, subgraphs)
mstgraph_lst = mstgraph_m$mstgraph # Updated MST 
subgraphs_m = mstgraph_m$subgraphs # Updated subgraphs
eid_btw_mst_m = mstgraph_m$eid_btw_mst # Edges between MSTs

csize_m = Rfast::Table(cluster_m) 
k_m = length(csize_m)

# Plot the updated MST
plotGraph(coords, mstgraph_lst) +
  geom_boundary(bnd_scaled) +
  labs(x = 'Coordinate X', y = 'Coordinate Y') +
  ggtitle('MST')

### Perform Change move
# first perform death move: (c1, c2) -> c2
merge_res = mergeCluster(mstgraph_lst, eid_btw_mst_m, subgraphs_m, csize_m, 
                         cluster_m, inc_mat, change = T)
# then perform birth move
split_res = splitCluster(mstgraph_lst, k_m-1, merge_res$subgraphs, merge_res$csize)

update_res_check_merge = updateMerge(merge_res, subgraphs_m, csize_m, eid_btw_mst_m, 
                                     cluster_m, edge_status, 
                                     adj_list, adj_edge_list, mstgraph_lst)

update_res_check_split = updateSplit(split_res, update_res_check_merge$subgraphs, k_m-1, update_res_check_merge$csize, 
                                     update_res_check_merge$eid_btw_mst, update_res_check_merge$cluster, 
                                     update_res_check_merge$estatus, adj_list, adj_edge_list)

cluster_change = update_res_check_split$cluster

# Create a temporary vector for the transformation (for the plot)
temp_cluster_change <- cluster_change
# Swap the assignments of 9 and 4
temp_cluster_change[cluster_change == 9] <- 8
temp_cluster_change[cluster_change == 8] <- 9
# Assign back the transformed values to cluster_change
cluster_change_toplot <- temp_cluster_change

plot_Graph_paper <- function(coords, graph, groups_assign_out, title = NULL) {
  require(ggplot2)
  
  # Extracting edge list from the graph
  edgelist <- get.edgelist(graph)
  
  # Preparing edge data
  edgedata <- data.frame(
    X1 = coords[edgelist[,1], 1],  # X coordinates of the starting points
    Y1 = coords[edgelist[,1], 2],  # Y coordinates of the starting points
    X2 = coords[edgelist[,2], 1],  # X coordinates of the ending points
    Y2 = coords[edgelist[,2], 2]   # Y coordinates of the ending points
  )
  
  # Preparing data for points: ensuring each unique coordinate is only plotted once
  pointsdata <- unique(rbind(
    data.frame(x = edgedata$X1, y = edgedata$Y1, Cluster = groups_assign_out[edgelist[, 1]]),
    data.frame(x = edgedata$X2, y = edgedata$Y2, Cluster = groups_assign_out[edgelist[, 2]])
  ))
  
  # Create the ggplot object with lines and points
  plot <- ggplot() +
    geom_segment(data = edgedata, aes(x = X1, y = Y1, xend = X2, yend = Y2), size = 1, color = "darkgrey") +
    geom_point(data = pointsdata, aes(x = x, y = y, color = as.factor(Cluster)), size = 4) +
    scale_color_brewer(palette = "Set1") +  # Using a Brewer color palette for distinct colors
    # labs(title = title, x = "X Coordinate", y = "Y Coordinate") +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_blank(),  # Remove y-axis title
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()
          # axis.ticks = element_line(color = "black"),
          # axis.line = element_line(colour = "black", size = 0.5,
          #                          arrow = arrow(type = "closed", length = unit(0.08, "inches")))
    )  # No y-axis title)
  return(plot)
}

# Define a color palette for spatial plots
my_palette_spat <- paletteer_d("ggthemes::Nuriel_Stone", 9)

# Plot the original and changed partitions
plot_change_old <- plot_Graph_paper(coords, mstgraph_lst, cluster_m ) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette_spat) +
  ggtitle(expression("Original Partition")) +
  theme(plot.title = element_text(size = 24))
plot_change_old

plot_change_new <- plot_Graph_paper(coords, mstgraph_lst, cluster_change_toplot ) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette_spat) +
  ggtitle(expression(bolditalic("Change Move"))) +
  theme(plot.title = element_text(size = 24))
plot_change_new

### Perform Hyper move
# update MST
mstgraph_new = proposeMST(graph0, edge_status, subgraphs_m)
mstgraph_lst_new = mstgraph_new$mstgraph

# Plot new MST
plot_new_mst <- plot_Graph_paper(coords, mstgraph_lst_new, cluster_m ) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette_spat) +
  ggtitle(expression(bolditalic("Hyper Move"))) +
  theme(plot.title = element_text(size = 24))
plot_new_mst


##### Save plots for the paper
combined_plots <- ggarrange(
  plot_change_new, plot_change_old,plot_new_mst, 
  ncol = 3, nrow = 1,
  labels = c("(a)", "(b)", "(c)"),
  label.x = -0.03, # Adjust the horizontal position of the labels
  label.y = 1.035, # Adjust the vertical position of the labels
  hjust = -0.6, # Horizontal justification of the labels
  vjust = 2 ,    # Vertical justification of the labels
  font.label = list(size = 22, face = "bold")
)
combined_plots

ggsave("./Plots/moves_explain_plots_paper.png", plot = combined_plots, 
       width = 20, height = 5, units = "in", dpi = 300, bg = "transparent")
