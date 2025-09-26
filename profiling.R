library(ggplot2)
library(igraph)

# Load Swiss-Roll Simulated Data
load("~/JASA supplementary/DP-RST/Simulations/Swiss_Roll_3k/Swiss_Roll_sim_3k_data_1.RData")

# set seed
set.seed(729)

# For debugging purposes
options(error=recover)

# Extract first three PCAs
Y_sample = sim_data[, 5:7]
loc = sim_data[, 1:2]
bnd = boundary

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

#Standartize the boundary given the params of coordinates
bnd_scaled = list()
bnd_scaled$x <- (bnd$x - mean(loc[,1]))/sd(loc[,1])
bnd_scaled$y <- (bnd$y - mean(loc[,2]))/sd(loc[,2])

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std)

##### SET UP INITIAL VALUES #####
# Get mesh and triangulation
mesh = gen2dMesh(coords, bnd_scaled)
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
mstgraph = mst(graph0)  # initial spanning tree
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph, 'weight')


# Define the "temperatures"
temp = c(1)
M = length(temp)

geodesic_dist_matrix <- distances(graph0, v = V(graph0), to = V(graph0),
                                  weights = E(graph0)$weight)


hc <- fastcluster::hclust(as.dist(geodesic_dist_matrix), method = "ward.D2")

# Cut the tree to form clusters
num_clusters <- 20
clusters <- cutree(hc, k = num_clusters)

# Compute centroids for each cluster in the batch
centroids <- aggregate(Y_std, by = list(clusters), FUN = mean)
centroids <- centroids[, -1]

# Initialize parameters for DP-RST

mu = list() # list for initial values of mu
sigmasq_y = list()
mstgraph_lst = list()  # initial spanning trees

clusters_bast <- clusters
clusters_bast_number <- length(unique(clusters_bast))
k_max = clusters_bast_number

cluster = matrix(rep(clusters_bast), nrow = n, ncol = M, byrow = F)
cluster_means_matrix <- as.matrix(centroids)
# # Get the means of each cluster
# for (i in 1:k_max) {
#   cluster_means_matrix[i,] <- colMeans(Y_std[c(which(clusters_bast == i)), , drop = FALSE])
# }

# Initialize all the groups from one team
teams = matrix(rep(1, k_max), ncol = M, nrow = k_max, byrow = F)
mu_teams = list()
# For each temperature, we have same initial cluster assignment
for(m in 1:M) {
  mu[[m]] = cluster_means_matrix #rep(0, p)
  mu_teams[[m]] = colMeans(Y_std)
  sigmasq_y[[m]] = cov(Y_std) #COvariance matrix
  mstgraph_lst[[m]] = mstgraph0
}

init_val = list()
init_val[['mstgraph_lst']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster
init_val[['teams']] = teams
init_val[['mu_teams']] = mu_teams
init_val[['sigmasq_y']] = sigmasq_y

hyperpar = list()
hyperpar[['sigmasq_mu']] = (0.5/(2*sqrt(1)))^2
hyperpar[['lambda_s']] = diag(1, p)
hyperpar[['nu']] = p
hyperpar[['M']] = M
hyperpar[['temp']] = temp
hyperpar[['k_max']] = k_max
hyperpar[['j_max']] = k_max

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 100 # MCMC iterations
BURNIN = 0 # burnin period length
THIN = 1     # thinning intervals

hyperpar[['alpha']] = 0.1 # Concentration parameter for the DPM model


library(profvis)
pv <- profvis({
  DP.RST(Y_std, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, PT = FALSE, seed = 6294)
})

# save for later comparison
htmlwidgets::saveWidget(pv, "profile_DP_RST.html", selfcontained = TRUE)


