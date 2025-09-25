##### LOAD NECESSARY LIBRARIES #####

library(DP.RST)
library(mclust)
library(ggplot2)
library(readr)
library(cluster)        # silhouette
library(clusterSim)     # Davies–Bouldin
library(clusterCrit)    # Calinski–Harabasz, Dunn

#-------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/DP-RST-workload")
getwd()

# Set seed
set.seed(7303)

# For debugging purposes
options(error=recover)
#-------------------------------------------------------------------------------
##### LOAD PRE-PROCESSED INTESTINE DATA #####

# Load pre-processed data
gut_df_wt_muscle <- readRDS("./Gut/gut_df_wt_muscle.rds")
load("./Gut/swiss_roll_wt_muscle_boundary.RData")

# Extract spatial coordinates and PCs
Y_sample = gut_df_wt_muscle[, 1:3]
loc = gut_df_wt_muscle[, c("x", "y")]
bnd = boundary

# Extract hystological labels
z = gut_df_wt_muscle$z

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

# Standardize the boundary to match the scaled coordinates
bnd_scaled = list()
bnd_scaled$x <- (bnd$x - mean(loc[,1]))/sd(loc[,1])
bnd_scaled$y <- (bnd$y - mean(loc[,2]))/sd(loc[,2])

# Create dataframe for plotting
df_subset <- data.frame(loc, z)

#-------------------------------------------------------------------------------
##### LOAD AND PROCESS DP-RST OUTPUT #####

load("./Gut/DP-RST_Outputs/Gut_DP.RST_FromNewBastPT_Version2_OutputOnly.RData")
output = Gut_DP.RST_FromNewBastPT_Version2

# Get the best partition
mode_based_partition = partition(DP.RST_output = output, method = "mode_based", 
                                 batch_size = 100)

#### Calculate the accuracy of the labeled observations #####
# Binary matrix for groups and teams
X <- table(sequence(length(mode_based_partition$groups_partition)), 
           mode_based_partition$groups_partition)
Z <- table(sequence(length(mode_based_partition$teams_partition)), 
           mode_based_partition$teams_partition)

# Get the membership of observations in each team
obs_in_teams <- X %*% Z
# Compute W such that w_ij = 1 if observation i and observation j share same team
obs_in_teams_vec <- obs_in_teams %*% sort(unique(mode_based_partition$teams_partition))

accur_teams = adjustedRandIndex(df_subset$z, obs_in_teams_vec) #
accur_teams # 0.4348815

DPM_partition <- obs_in_teams_vec

#-------------------------------------------------------------------------------
##### LOAD BAYES-SPACE RESULTS #####

BayesSpace_Gut_clustering_results <- readRDS("./Gut/BayesSpace_Gut_clustering_results.rds")

BS.df <- data.frame(x = BayesSpace_Gut_clustering_results$imagerow, 
                    y = BayesSpace_Gut_clustering_results$imagecol, 
                    BS_z = BayesSpace_Gut_clustering_results$cluster_q5_pc3)

BS_merged_df <- merge(BS.df, df_subset,
                      by = c("x", "y"),
                      all = TRUE)

BS_accur = adjustedRandIndex(BS_merged_df$z, BS_merged_df$BS_z)
BS_accur # 0.1280336

# Reorder the partition based on a reference vector
BS_partition <- BS_merged_df$BS_z

#-------------------------------------------------------------------------------
##### LOAD SC-MEB RESULTS #####

SC_MEB_Gut_clustering_results <- readRDS("./Gut/SC-MEB_Gut_clustering_results.rds")

SC_MEB.df <- data.frame(x = SC_MEB_Gut_clustering_results$x, 
                        y = SC_MEB_Gut_clustering_results$y, 
                        SC_z = SC_MEB_Gut_clustering_results$pc3_k5)

SC_MEB_merged_df <- merge(SC_MEB.df, df_subset,
                          by = c("x", "y"),
                          all = TRUE)

SC.MEB_accur = adjustedRandIndex(SC_MEB_merged_df$z, SC_MEB_merged_df$SC_z)
SC.MEB_accur # 0.1124159

# Reorder the partition based on a reference vector
SC.MEB_partition <- SC_MEB_merged_df$SC_z

#-------------------------------------------------------------------------------
##### LOAD DR-SC RESULTS #####

DR_SC_Gut_clustering_results <- readRDS("./Gut/DR-SC_Gut_clustering_results.rds")

DR_SC.df <- data.frame(x = DR_SC_Gut_clustering_results$imagerow, 
                       y = DR_SC_Gut_clustering_results$imagecol, 
                       DR_z = DR_SC_Gut_clustering_results$cluster_q5_pc3)

DR_SC_merged_df <- merge(DR_SC.df, df_subset,
                         by = c("x", "y"),
                         all = TRUE)

DR_SC_accur = adjustedRandIndex(DR_SC_merged_df$z, DR_SC_merged_df$DR_z)
DR_SC_accur # 0.1530579

# Reorder the partition based on a reference vector
DR.SC_partition <- DR_SC_merged_df$DR_z

#-------------------------------------------------------------------------------
##### LOAD K-MEANS RESULTS #####

kmeans_Gut_clustering_results <- readRDS("./Gut/kmeans_Gut_clustering_results.rds")

kmeans.df <- data.frame(x = kmeans_Gut_clustering_results$x, 
                        y = kmeans_Gut_clustering_results$y, 
                        kmeans_z = kmeans_Gut_clustering_results$pc3_k5)

kmeans_merged_df <- merge(kmeans.df, df_subset,
                          by = c("x", "y"),
                          all = TRUE)

kmeans_accur = adjustedRandIndex(kmeans_merged_df$z, kmeans_merged_df$kmeans_z)
kmeans_accur # 0.02973941

# Reorder the partition based on a reference vector
kmeans_partition <- kmeans_merged_df$kmeans_z

#-------------------------------------------------------------------------------
##### LOAD SEDR RESULTS #####

sedr_gut_results <- read_csv("./Competing_Methods/SEDR_runs/res_realdata_csv/sedr_gut_reduced_pca3_results.csv")

sedr.df <- data.frame(y = sedr_gut_results$x, 
                      x = sedr_gut_results$y, 
                      sedr_z = sedr_gut_results$mclust_pca3_5)

sedr_merged_df <- merge(sedr.df, df_subset,
                        by = c("x", "y"),
                        all = TRUE)

sedr_accur = adjustedRandIndex(sedr_merged_df$z, sedr_merged_df$sedr_z)
sedr_accur # 0.0481976

# Reorder the partition based on a reference vector
sedr_partition <- sedr_merged_df$sedr_z

#-------------------------------------------------------------------------------
##### LOAD GRAPHST RESULTS #####

graphst_gut_results <- read_csv("./Competing_Methods/GraphST_runs/res_realdata_csv/graphst_gut_reduced_results.csv")

graphst.df <- data.frame(y = graphst_gut_results$x, 
                         x = graphst_gut_results$y, 
                         graphst_z = graphst_gut_results$mclust_5)

graphst_merged_df <- merge(graphst.df, df_subset,
                           by = c("x", "y"),
                           all = TRUE)

graphst_accur = adjustedRandIndex(graphst_merged_df$z, graphst_merged_df$graphst_z)
graphst_accur # 0.1419926

# Reorder the partition based on a reference vector
graphst_partition <- graphst_merged_df$graphst_z

#-------------------------------------------------------------------------------
##### LOAD STAGATE RESULTS #####

stagate_gut_results <- read_csv("./Competing_Methods/STAGATE_runs/res_realdata_csv/stagate_gut_reduced_results.csv")

stagate.df <- data.frame(y = stagate_gut_results$x, 
                         x = stagate_gut_results$y, 
                         stagate_z = stagate_gut_results$mclust_5)

stagate_merged_df <- merge(stagate.df, df_subset,
                           by = c("x", "y"),
                           all = TRUE)

stagate_accur = adjustedRandIndex(stagate_merged_df$z, stagate_merged_df$stagate_z)
stagate_accur # 0.2669045

# Reorder the partition based on a reference vector
stagate_partition <- stagate_merged_df$stagate_z

#-------------------------------------------------------------------------------
##### LOAD SpaGCN RESULTS #####

spagcn_gut_results <- read_csv("Competing_Methods/SpaGCN_runs/res_realdata_csv/all_clustering_results_Gut_reduced.csv")

spagcn.df <- data.frame(spots = spagcn_gut_results$...1, 
                        x = spagcn_gut_results$array_row, 
                        y = spagcn_gut_results$array_col, 
                        spagcn_z_with = (spagcn_gut_results$refined_pred_5clusters_3pcs_with_histology+1),
                        spagcn_z_without = (spagcn_gut_results$refined_pred_5clusters_3pcs_without_histology+1))

df_subset$spots <- rownames(df_subset)
spagcn_merged_df <- merge(df_subset, spagcn.df, by = "spots")

spagcn_accur_with = adjustedRandIndex(spagcn_merged_df$spagcn_z_with, spagcn_merged_df$z)
spagcn_accur_with # 0.1304065

# Reorder the partition based on a reference vector
spagcn_partition_with <- spagcn_merged_df$spagcn_z_with

spagcn_accur_without = adjustedRandIndex(spagcn_merged_df$spagcn_z_without, spagcn_merged_df$z)
spagcn_accur_without # 0.1990036

# Reorder the partition based on a reference vector
spagcn_partition_without <- spagcn_merged_df$spagcn_z_without

#-------------------------------------------------------------------------------
##### LOAD stLearn RESULTS #####

stlearn_gut_results_with <- read_csv("./Competing_Methods/stLearn_runs/res_realdata_csv/Gut_reduced_clustering_3PCs_5clusters_with_image.csv")
stlearn_gut_results_without <- read_csv("./Competing_Methods/stLearn_runs/res_realdata_csv/Gut_reduced_clustering_3PCs_5clusters_no_image.csv")

stlearn.df <- data.frame(y = stlearn_gut_results_with$x, 
                         x = stlearn_gut_results_with$y, 
                         stlearn_z_with = (stlearn_gut_results_with$cluster_5+1),
                         stlearn_z_without = (stlearn_gut_results_without$cluster_5+1))

stlearn_merged_df <- merge(stlearn.df, df_subset,
                           by = c("x", "y"),
                           all = TRUE)

stlearn_accur_with = adjustedRandIndex(stlearn_merged_df$stlearn_z_with, stlearn_merged_df$z)
stlearn_accur_with # 0.2753501

# Reorder the partition based on a reference vector
stlearn_partition_with <- stlearn_merged_df$stlearn_z_with

stlearn_accur_without = adjustedRandIndex(stlearn_merged_df$stlearn_z_without, stlearn_merged_df$z)
stlearn_accur_without # 0.1507263

# Reorder the partition based on a reference vector
stlearn_partition_without <- stlearn_merged_df$stlearn_z_without

#-------------------------------------------------------------------------------

##### 1.  Helper to collect partitions ########################################
partitions <- list(
  DP_RST            = DPM_partition,
  BayesSpace        = BS_partition,
  SC_MEB            = SC.MEB_partition,
  DR_SC             = DR.SC_partition,
  k_means           = kmeans_partition,
  SEDR              = sedr_partition,
  GraphST           = graphst_partition,
  STAGATE           = stagate_partition,
  SpaGCN_with_img   = spagcn_partition_with,
  SpaGCN_no_img     = spagcn_partition_without,
  stLearn_with_img  = stlearn_partition_with,
  stLearn_no_img    = stlearn_partition_without
)

##### 2.  Metric engine #######################################################
compute_metrics <- function(X, cl) {
  cl            <- as.integer(factor(cl))     # ensure 1…K
  n             <- nrow(X)
  k             <- length(unique(cl))
  
  ## Within-cluster variation
  centroids     <- rowsum(X, cl) / as.vector(table(cl))
  WSS           <- sum(rowSums((X - centroids[cl, ])^2))
  avg_within    <- WSS / n
  
  ## Between-cluster variation
  if (k > 1) {
    d_cent       <- as.matrix(dist(centroids))^2
    avg_between  <- mean(d_cent[upper.tri(d_cent)])
  } else {
    avg_between  <- NA
  }
  
  ## Objective indices --------------------------------------------------------
  if (k > 1) {
    sil          <- mean(cluster::silhouette(cl, dist(X))[, "sil_width"])
    dbi          <- clusterSim::index.DB(X, cl, centrotypes = "centroids")$DB
    chi          <- clusterCrit::intCriteria(X, cl, "Calinski_Harabasz")$calinski_harabasz
    dunn         <- clusterCrit::intCriteria(X, cl, "Dunn")$dunn
  } else {
    sil <- dbi <- chi <- dunn <- NA
  }
  
  c(avg_within = avg_within,
    avg_between = avg_between,
    silhouette = sil,
    davies_bouldin = dbi,
    calinski_harabasz = chi,
    dunn = dunn)
}




##### 3.  Assemble results ####################################################
Y_mat     <- as.matrix(Y_std)            # 3-column (PC1–PC3) matrix
metrics_df <- do.call(rbind, lapply(names(partitions), function(meth) {
  res <- compute_metrics(Y_mat, partitions[[meth]])
  data.frame(Method = meth, t(res), row.names = NULL, check.names = FALSE)
}))
print(metrics_df, digits = 3)

# Optional: write.csv(metrics_df, "cluster_evaluation_metrics.csv", row.names = FALSE)








