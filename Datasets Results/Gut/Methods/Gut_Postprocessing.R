library(ggplot2)
library(igraph)
library(ggpubr)
library(patchwork)

#-------------------------------------------------------------------------------

load("Gut/swiss_roll_wt_muscle_finaltouches1.RData")
load("Gut/swiss_roll_wt_muscle_boundary.RData")

# Extract first three PCAs
Y_sample = swiss_roll_wt_muscle_finaltouches1[, 1:3]
loc = swiss_roll_wt_muscle_finaltouches1[, 4:5]

bnd = boundary

z = swiss_roll_wt_muscle_finaltouches1$z
z_man = swiss_roll_wt_muscle_finaltouches1$z_man

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

#Standartize the boundary given the params of coordinates
bnd_scaled = list()
bnd_scaled$x <- (bnd$x - mean(loc[,1]))/sd(loc[,1])
bnd_scaled$y <- (bnd$y - mean(loc[,2]))/sd(loc[,2])

#-------------------------------------------------------------------------------

file_path <- "Gut/DP-RST_Outputs/Gut_DP.RST_FromNewBastPT_p10_Version2_OutputOnly.RData"
load(file_path)
output = Gut_DP.RST_FromNewBastPT_p10_Version2


file_path_bast <- "Gut/DP-RST_Outputs/Gut_BAST_PT_OutputOnly.RData"
load(file_path_bast)
output_bast = Gut_BAST_PT

#-------------------------------------------------------------------------------

marg_lik <- ggplot(data.frame(x = seq_along(output[["marginal_likelihood_out"]]),
                              y = output[["marginal_likelihood_out"]]), aes(x, y)) +
  geom_line() +
  ggtitle("Marginal Likelihood post BurnIn and Thinning") +
  xlab("Iteration") +
  ylab("values")
marg_lik

#-------------------------------------------------------------------------------

table(output[["moves_track_out"]])

table(output[["j_teams_out"]])
plot(output[["j_teams_out"]])

teams_freq <- as.data.frame(table(output[["j_teams_out"]]))

moves_freq_plot <- ggplot(teams_freq, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "grey") +
  labs(title = "Frequency of # of refined clusters", x = "# of Teams", y = "Frequency\n") +
  theme_classic()
moves_freq_plot

#-------------------------------------------------------------------------------

accur_spatial <- c()
accur_refined <- c()

for (i in 1:length(output[["marginal_likelihood_out"]])){
  Gut_df <- data.frame(res_memb = output[["cluster_out"]][i,], z)
  Gut_df <- subset(Gut_df, !is.na(z))

  accur_spatial[i] <- mclust::adjustedRandIndex(Gut_df$res_memb,
                                                Gut_df$z)

  X <- table(sequence(length(output[["cluster_out"]][i,])),
             output[["cluster_out"]][i,])
  Z <- table(sequence(length(output[["teams_out"]][[i]])),
             output[["teams_out"]][[i]])

  # Get the membership of observations in each team
  obs_in_teams <- X %*% Z
  # Compute W such that w_ij = 1 if observation i and observation j share same team
  obs_in_teams_vec <- obs_in_teams %*% sort(unique(output[["teams_out"]][[i]]))

  Gut_df <- data.frame(res_memb = obs_in_teams_vec, z)
  Gut_df <- subset(Gut_df, !is.na(z))

  accur_refined[i] = mclust::adjustedRandIndex(Gut_df$res_memb, Gut_df$z)
}

acc_spatial <- ggplot(data.frame(x = seq_along(accur_spatial), y = accur_spatial), aes(x, y)) +
  geom_line() +
  ggtitle("Accuracy of Spatial Partition") +
  xlab("Iteration") +
  ylab("Adjusted Rand Index (ARI)")
acc_spatial

acc_refined <- ggplot(data.frame(x = seq_along(accur_refined), y = accur_refined), aes(x, y)) +
  geom_line() +
  ggtitle("Accuracy of Refined Partition") +
  xlab("Iteration") +
  ylab("Adjusted Rand Index (ARI)")
acc_refined

#-------------------------------------------------------------------------------

# Get the best partition
mode_based_partition = partition(DP.RST_output = output, method = "mode_based", batch_size = 100)
frobenius_partition = partition(DP.RST_output = output, method = "frobenius", batch_size = 100)

#-------------------------------------------------------------------------------

### Accuracy of the crude spatial partition
Gut_df <- data.frame(res_1 = mode_based_partition$groups_partition,
                     res_2 = frobenius_partition$groups_partition,
                     z)
Gut_df <- subset(Gut_df, !is.na(z))



## Accuracy of the refined partition for the mode_based method
# Binary membership matrix for groups and teams
X <- table(sequence(length(mode_based_partition$groups_partition)),
           mode_based_partition$groups_partition)
Z <- table(sequence(length(mode_based_partition$teams_partition)),
           mode_based_partition$teams_partition)

# Get the membership of observations in each team
obs_in_teams <- X %*% Z
# Compute W such that w_ij = 1 if observation i and observation j share same team
obs_in_teams_vec_mode_based <- obs_in_teams %*% sort(unique(mode_based_partition$teams_partition))



## Accuracy of the refined partition for the mode_based method
# Binary membership matrix for groups and teams
X <- table(sequence(length(frobenius_partition$groups_partition)),
           frobenius_partition$groups_partition)
Z <- table(sequence(length(frobenius_partition$teams_partition)),
           frobenius_partition$teams_partition)

# Get the membership of observations in each team
obs_in_teams <- X %*% Z
# Compute W such that w_ij = 1 if observation i and observation j share same team
obs_in_teams_vec_frobenius <- obs_in_teams %*% sort(unique(frobenius_partition$teams_partition))


Gut_df <- data.frame(res_groups_1 = mode_based_partition$groups_partition,
                     res_teams_1 = obs_in_teams_vec_mode_based,
                     res_groups_2 = frobenius_partition$groups_partition,
                     res_teams_2 = obs_in_teams_vec_frobenius,
                     z, z_man)
Gut_df <- subset(Gut_df, !is.na(z))


accur_groups_1 = mclust::adjustedRandIndex(Gut_df$res_groups_1, Gut_df$z)
accur_groups_1

accur_groups_2 = mclust::adjustedRandIndex(Gut_df$res_groups_2, Gut_df$z)
accur_groups_2

accur_teams_1 = mclust::adjustedRandIndex(Gut_df$res_teams_1, Gut_df$z)
accur_teams_1

accur_teams_2 = mclust::adjustedRandIndex(Gut_df$res_teams_2, Gut_df$z)
accur_teams_2

#-------------------------------------------------------------------------------

### Plot crude spatial partition

groups_plotted_mode_based <- groups_plot(coords, mode_based_partition$groups_partition, bnd_scaled) +
  theme(legend.position = "none") +
  ggtitle(paste("DP-RST Spatial Partition (mode), k = ", length(unique(mode_based_partition$groups_partition)), " , ARI = ",
                round(accur_groups_1, 3),
                sep = ""))


groups_plotted_frobenius <- groups_plot(coords, frobenius_partition$groups_partition, bnd_scaled) +
  theme(legend.position = "none") +
  ggtitle(paste("DP-RST Spatial Partition (frobenius), k = ", length(unique(frobenius_partition$groups_partition)), " , ARI = ",
                round(accur_groups_2, 3),
                sep = ""))

#-------------------------------------------------------------------------------

### Plot the refined partition

teams_plotted_mode_based <- teams_plot(coords, mode_based_partition$teams_partition,
                            mode_based_partition$groups_partition, bnd_scaled) +
  theme(legend.position = "none") +
  ggtitle(paste("DP-RST Refined Partition (mode), k = ", length(unique(mode_based_partition$teams_partition)), " , ARI = ",
                round(accur_teams_1, 3),
                sep = ""))


teams_plotted_frobenius <- teams_plot(coords, frobenius_partition$teams_partition,
                                       frobenius_partition$groups_partition, bnd_scaled) +
  theme(legend.position = "none") +
  ggtitle(paste("DP-RST Refined Partition (frobenius), k = ", length(unique(frobenius_partition$teams_partition)), " , ARI = ",
                round(accur_teams_2, 3),
                sep = ""))

#-------------------------------------------------------------------------------

true_clusters <- groups_plot(coords, as.factor(z), bnd_scaled) +
  theme(legend.position = "none") +
  ggtitle(paste("True cluster memberships, k = ", length(unique(z[!is.na(z)])), sep = ""))
# true_clusters

#-------------------------------------------------------------------------------
### Get the initial partition supplied to DP-RST

##### CHOOSE THE ITERATION #####
groups_assign = output_bast$cluster_out; groups_number = output_bast$k_out
# List to save co-clustering matrices
# List to save co-clustering matrices
W_cum <- matrix(0, nrow = n, ncol = n)

# Get the mode of the team
tbl <- table(groups_number)
max_freq <- max(tbl)
modes <- as.numeric(names(tbl[tbl == max_freq]))

# Get the indices where teams_number equals the mode
mode_indices <- which(groups_number == modes)

# Run the loop to create co-clustering matrix for over the selected indices
for (i in 1:length(mode_indices)) {
  # Binary matrix for groups and teams
  X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
  # Compute W such that w_ij = 1 if observation i and observation j share same team
  W <- X %*% t(X)
  W_cum <- W_cum + W
}

mean_matrix <- W_cum/length(mode_indices)

# Calculate the vector with Frobenius norms
norm = c()

for (i in 1:length(mode_indices)) {
  # Binary matrix for groups and teams
  X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
  # Get the membership of observations in each team
  mat <- X %*% t(X)
  diff_matrix <- mat - mean_matrix
  norm[i] <- norm(diff_matrix, type = "F")
}

# Fins the minimum Frobenius norm index
min_index <- mode_indices[which.min(norm)]
groups_assign_out <- groups_assign[min_index, ]
rm(W, W_cum, mean_matrix)

#-------------------------------------------------------------------------------

groups_plotted_BASTinit <- groups_plot(coords, groups_assign_out, bnd_scaled) +
  theme(legend.position = "none") +
  ggtitle(paste("Initialization from BAST PT, k = ", length(unique(groups_assign_out)),
                sep = ""))

#-------------------------------------------------------------------------------

# Extract hyperparameters dynamically
sigmasq_mu <- output$hyperpar$sigmasq_mu
nu <- output$hyperpar$nu
alpha <- output$hyperpar$alpha
M <- output$hyperpar$M
k_max <- output$hyperpar$k_max
j_max <- output$hyperpar$j_max
PCs <- length(output$init_val$mu_teams[[1]])

#-------------------------------------------------------------------------------

# Extract the file name dynamically (removes path and extension)
file_name <- basename(file_path)
file_name <- sub("\\.RData$", "", file_name)  # Remove the .RData extension

file_name_bast <- basename(file_path_bast)
file_name_bast <- sub("\\.RData$", "", file_name_bast)  # Remove the .RData extension

#-------------------------------------------------------------------------------

# Create the dynamic title
title_text <- sprintf("Gut, DP-RST with PT, PCs = %d, diff = 0.1. Initialize at BAST PT (new initialization).
sigmasq_mu = %.3f, nu = %d, alpha = %.2f, M = %d, k_max = %d, j_max = %d. MCMC = 25,000; BURNIN = 20,000; THIN = 1.
Model seed 5284. File name %s . File name of initialization %s .", PCs, sigmasq_mu, nu, alpha, M, k_max, j_max, file_name, file_name_bast)

#-------------------------------------------------------------------------------

# Update the plot annotation
plots <- ggarrange(true_clusters, teams_plotted_mode_based, teams_plotted_frobenius,
                   groups_plotted_BASTinit, groups_plotted_mode_based, groups_plotted_frobenius,
                   moves_freq_plot, acc_refined, marg_lik,
                   ncol = 3, nrow = 3) +
  plot_annotation(title = title_text,
                  theme = theme(plot.title = element_text(size = 14, color="darkblue", face = 'bold')))

#-------------------------------------------------------------------------------

# Extract the base file name and remove "OutputOnly"
save_file_name <- sub("_OutputOnly", "", file_name)  # Removes "OutputOnly" if present
save_path <- paste0("Gut/DP-RST_Outputs/", save_file_name, ".pdf")

ggsave(file = save_path, plots, width = 8, height = 6, scale = 2)

#-------------------------------------------------------------------------------
