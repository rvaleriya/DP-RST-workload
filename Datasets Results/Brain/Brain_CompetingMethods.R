library(ggplot2)
library(dplyr)
library(purrr)
library(mclust)
library(tidyr)
library(stringr)
#-------------------------------------------------------------------------------
set.seed(42)
### BayesSpace ###
library(BayesSpace)
sce_brain <- readVisium("Brain/Space_Ranger_Data_Brain")

brain <- spatialPreprocess(sce_brain)

# brain_tune <- qTune(brain, qs=seq(2, 20), d=3)
# qPlot(brain_tune)

library(BayesSpace)

clusters <- c(5, 7, 9)

# Preallocate lists to store results
results_q3 <- list()
results_q10 <- list()

for (c in clusters) {
  # Clustering with q=3
  brain_q3 <- spatialCluster(brain, q=c, d=3, platform='Visium',
                             init.method="mclust", model="t", gamma=3,
                             nrep=50000, burn.in=40000, save.chain=FALSE)

  # Store clustering labels
  df_q3 <- data.frame(
    imagerow = brain_q3$imagerow,
    imagecol = brain_q3$imagecol,
    spatial_cluster = brain_q3$spatial.cluster,
    PCs = 3,
    q = c
  )
  results_q3[[as.character(c)]] <- df_q3

  # Clustering with q=10
  brain_q10 <- spatialCluster(brain, q=c, d=10, platform='Visium',
                              init.method="mclust", model="t", gamma=3,
                              nrep=50000, burn.in=40000, save.chain=FALSE)

  # Store clustering labels
  df_q10 <- data.frame(
    imagerow = brain_q10$imagerow,
    imagecol = brain_q10$imagecol,
    spatial_cluster = brain_q10$spatial.cluster,
    PCs = 10,
    q = c
  )
  results_q10[[as.character(c)]] <- df_q10
}


rename_with_id <- function(df, q, pc) {
  colname <- paste0("cluster_q", q, "_pc", pc)
  df %>%
    select(imagerow, imagecol, spatial_cluster) %>%
    rename_with(~ colname, .cols = spatial_cluster)
}

# Generate merged results for PC=3
df_list_q3 <- lapply(names(results_q3), function(q) {
  rename_with_id(results_q3[[q]], q = q, pc = 3)
})

# Generate merged results for PC=10
df_list_q10 <- lapply(names(results_q10), function(q) {
  rename_with_id(results_q10[[q]], q = q, pc = 10)
})

# Merge all dataframes by coordinates
merged_df <- purrr::reduce(c(df_list_q3, df_list_q10), full_join, by = c("imagerow", "imagecol"))
saveRDS(merged_df, file = "Brain/BayesSpace_Brain_clustering_results.rds")

load("Brain/DFPLC_151510_data_for_iIMPACT.RData")

merged_df_renamed <- merged_df %>%
  rename(x = imagerow, y = imagecol)

brain_merged_df <- merge(ground_truth_df, merged_df_renamed,
                         by = c("x", "y"),
                         all = TRUE)

# brain_BS_labels <- ggplot() +
#   # geom_boundary(bnd_scaled_r) +
#   geom_point(aes(x = x, y = y, col = as.factor(z_BS)), data = brain_merged_df) +
#   labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "Clusters") +
#   ggtitle('BayesSpace Labels') +
#   theme(legend.position = "bottom")
# brain_BS_labels

brain_subset <- subset(brain_merged_df, !is.na(brain_merged_df$type))
head(brain_subset)

clustering_cols <- grep("^cluster_", names(brain_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  adjustedRandIndex(brain_subset$type, brain_subset[[col]])
})

# Extract q and PC from clustering column names
ari_df <- data.frame(Clustering = names(ari_scores), ARI = ari_scores) %>%
  mutate(
    q = as.integer(str_extract(Clustering, "(?<=q)\\d+")),
    PC = as.integer(str_extract(Clustering, "(?<=pc)\\d+"))
  ) %>%
  arrange(PC, q)

# Split into two tables
ari_pc3 <- ari_df %>% filter(PC == 3) %>% arrange(desc(ARI))
ari_pc10 <- ari_df %>% filter(PC == 10) %>% arrange(desc(ARI))

# View results
print(ari_pc3)

# Clustering       ARI q PC
# cluster_q7_pc3 cluster_q7_pc3 0.3559351 7  3
# cluster_q9_pc3 cluster_q9_pc3 0.3200856 9  3
# cluster_q5_pc3 cluster_q5_pc3 0.2870042 5  3

print(ari_pc10)

# Clustering       ARI q PC
# cluster_q5_pc10 cluster_q5_pc10 0.5037085 5 10
# cluster_q7_pc10 cluster_q7_pc10 0.4196511 7 10
# cluster_q9_pc10 cluster_q9_pc10 0.3831041 9 10

#-------------------------------------------------------------------------------
##### SC-MEB #####
set.seed(42)

library(SC.MEB)
load("Brain/DFPLC_151510_data_for_iIMPACT.RData")

### Set the parameters
#platform = "ST"
beta_grid = seq(0, 5, 0.2) #Vector specifying the smoothness of Random Markov Field
K_set= c(5, 7, 9) # 2:10 #Numbers of mixture components
parallel=TRUE #Logical value specifing the run the model in parallel
num_core = 3
PX = TRUE #Logical value for paramter expansion in EM algorithm
maxIter_ICM = 10 #Maximum iteration of ICM algorithm
maxIter = 50 #Maximum iteration of EM algorithm

#### Calculating the neighborhood from coordinates
Adj_sp <- getneighborhood_fast(as.matrix(loc), cutoff = 1.2)

#### Run the SC-MEB in parallel
brain_SC.MEB_3pc = SC.MEB(Y[, 1:3], Adj_sp, beta_grid = beta_grid, K_set = K_set,
                    parallel = parallel, num_core = num_core, PX = PX,
                    maxIter_ICM = maxIter_ICM, maxIter = maxIter)
str(brain_SC.MEB_3pc[,1])
str(brain_SC.MEB_3pc[,2])
str(brain_SC.MEB_3pc[,3])

brain_SC.MEB_10pc = SC.MEB(Y[, 1:10], Adj_sp, beta_grid = beta_grid, K_set = K_set,
                          parallel = parallel, num_core = num_core, PX = PX,
                          maxIter_ICM = maxIter_ICM, maxIter = maxIter)

# The item 'x' is clustering label.
# The item 'ell' is the opposite log-likelihood for each beta and K.
# The item 'mu' is the mean of each component.
# The item 'sigma' is the variance of each component.
# The item 'gam' is the posterior probability.
# The item 'beta' is the estimated smoothing parameter.


# #### Selecting the number of clusters using BIC
# selectKPlot(brain_fit, K_set = K_set, criterion = "BIC")
#
# #### Selecting the number of clusters using Modified BIC
# selectKPlot(brain_fit, K_set = K_set, criterion = "MBIC")
#
# out = selectK(brain_fit, K_set = K_set, criterion = "BIC")
# ClusterPlot(out, loc)
#
# brain_SC.MEB_df <- data.frame(loc, out[["best_K_label"]], z)
z = ground_truth_df$type

brain_SC.MEB_df <- data.frame(loc,
                              pc3_k5 = brain_SC.MEB_3pc[,1]$x,
                              pc3_k7 = brain_SC.MEB_3pc[,2]$x,
                              pc3_k9 = brain_SC.MEB_3pc[,3]$x,
                              pc10_k5 = brain_SC.MEB_10pc[,1]$x,
                              pc10_k7 = brain_SC.MEB_10pc[,2]$x,
                              pc10_k9 = brain_SC.MEB_10pc[,3]$x,
                              z)

saveRDS(brain_SC.MEB_df, file = "Brain/SC-MEB_Brain_clustering_results.rds")

brain_subset <- subset(brain_SC.MEB_df, !is.na(z))
head(brain_subset)

clustering_cols <- grep("^pc", names(brain_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  round(adjustedRandIndex(brain_subset$z, brain_subset[[col]]), 3)
})

# Make dataframe and extract PC + K info
ari_df <- data.frame(
  Clustering = names(ari_scores),
  ARI = ari_scores
) %>%
  mutate(
    PC = str_extract(Clustering, "pc\\d+"),
    K = as.integer(str_extract(Clustering, "(?<=k)\\d+"))
  ) %>%
  arrange(PC, K)

# Split into two tables
ari_pc3 <- ari_df %>% filter(PC == "pc3") %>% arrange(desc(ARI))
ari_pc10 <- ari_df %>% filter(PC == "pc10") %>% arrange(desc(ARI))

# View results
print(ari_pc3)

# Clustering   ARI  PC K
# pc3_k5     pc3_k5 0.358 pc3 5
# pc3_k7     pc3_k7 0.263 pc3 7
# pc3_k9     pc3_k9 0.238 pc3 9

print(ari_pc10)
# Clustering   ARI   PC K
# pc10_k5    pc10_k5 0.346 pc10 5
# pc10_k7    pc10_k7 0.264 pc10 7
# pc10_k9    pc10_k9 0.227 pc10 9


ggplot(brain_subset, aes(x = x, y = y, color = as.factor(pc3_k7))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=7")



#--------------------------------------------------------------------------------
##### DR-SC #####
set.seed(42)

library(Seurat)
library(DR.SC)
### Create Seurat object ###
# Create count matrix
data_dir_brain_count <- 'Brain/Space_Ranger_Data_brain/filtered_feature_bc_matrix'
list.files(data_dir_brain_count) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir_brain_count)

# Load spatial data
tissue_positions_list <- readr::read_csv("Brain/Space_Ranger_Data_brain/spatial/tissue_positions_list.csv",
                                  col_names = FALSE)
head(tissue_positions_list)
dim(tissue_positions_list)

meta_data <- data.frame(tissue_positions_list[, c(1,5,6)])
# Keep only rows in meta_data where X1 values are in the column names of expression_matrix
meta_data <- meta_data[meta_data$X1 %in% colnames(expression_matrix), ]
head(meta_data)
# Reorder rows of meta_data based on the order of X1 in colnames(expression_matrix)
meta_data <- meta_data[match(colnames(expression_matrix), meta_data$X1), ]
head(meta_data)

# Rename row names as X1
rownames(meta_data) <- meta_data$X1
# Delete the first column
meta_data <- meta_data[, -1]
# Rename the columns of meta_data
colnames(meta_data) <- c("row", "col")

# row.names(meta_data) <- colnames(expression_matrix)
head(meta_data)

# Create Seurat object
seurat_object = CreateSeuratObject(counts = expression_matrix, meta.data = meta_data)
head(seurat_object)

# standard log-normalization
brain <- NormalizeData(seurat_object, verbose = F)
# choose 2000 highly variable features
seu <- FindVariableFeatures(brain, nfeatures = 2000, verbose = F)

# Fit DR-SC model using 2000 highly variable features
clusters <- c(5, 7, 9)

# Preallocate lists to store results
results_q3 <- list()
results_q10 <- list()

for (c in clusters) {
  seu_fit_pc3 <- DR.SC(seu, K=c, q = 3, platform = 'Visium', verbose=F)

  # Store clustering labels
  df_q3 <- data.frame(
    imagerow = seu_fit_pc3$row,
    imagecol = seu_fit_pc3$col,
    spatial_cluster = seu_fit_pc3$spatial.drsc.cluster,
    PCs = 3,
    q = c
  )
  results_q3[[as.character(c)]] <- df_q3

  seu_fit_pc10 <- DR.SC(seu, K=c, q = 10, platform = 'Visium', verbose=F)

  df_q10 <- data.frame(
    imagerow = seu_fit_pc10$row,
    imagecol = seu_fit_pc10$col,
    spatial_cluster = seu_fit_pc10$spatial.drsc.cluster,
    PCs = 3,
    q = c
  )
  results_q10[[as.character(c)]] <- df_q10
}


rename_with_id <- function(df, q, pc) {
  colname <- paste0("cluster_q", q, "_pc", pc)
  df %>%
    select(imagerow, imagecol, spatial_cluster) %>%
    rename_with(~ colname, .cols = spatial_cluster)
}

# Generate merged results for PC=3
df_list_q3 <- lapply(names(results_q3), function(q) {
  rename_with_id(results_q3[[q]], q = q, pc = 3)
})

# Generate merged results for PC=10
df_list_q10 <- lapply(names(results_q10), function(q) {
  rename_with_id(results_q10[[q]], q = q, pc = 10)
})

# Merge all dataframes by coordinates
merged_df <- purrr::reduce(c(df_list_q3, df_list_q10), full_join, by = c("imagerow", "imagecol"))
saveRDS(merged_df, file = "Brain/DR-SC_Brain_clustering_results.rds")


load("Brain/DFPLC_151510_data_for_iIMPACT.RData")

merged_df_renamed <- merged_df %>%
  rename(x = imagerow, y = imagecol)

brain_merged_df <- merge(ground_truth_df, merged_df_renamed,
                         by = c("x", "y"),
                         all = TRUE)


brain_subset <- subset(brain_merged_df, !is.na(brain_merged_df$type))
head(brain_subset)

clustering_cols <- grep("^cluster_", names(brain_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  adjustedRandIndex(brain_subset$type, brain_subset[[col]])
})

# Extract q and PC from clustering column names
ari_df <- data.frame(Clustering = names(ari_scores), ARI = ari_scores) %>%
  mutate(
    q = as.integer(str_extract(Clustering, "(?<=q)\\d+")),
    PC = as.integer(str_extract(Clustering, "(?<=pc)\\d+"))
  ) %>%
  arrange(PC, q)

# Split into two tables
ari_pc3 <- ari_df %>% filter(PC == 3) %>% arrange(desc(ARI))
ari_pc10 <- ari_df %>% filter(PC == 10) %>% arrange(desc(ARI))

# View results
print(ari_pc3)
# Clustering       ARI q PC
# cluster_q5_pc3 cluster_q5_pc3 0.2458759 5  3
# cluster_q9_pc3 cluster_q9_pc3 0.2395309 9  3
# cluster_q7_pc3 cluster_q7_pc3 0.2366944 7  3

print(ari_pc10)
# Clustering       ARI q PC
# cluster_q5_pc10 cluster_q5_pc10 0.3176150 5 10
# cluster_q7_pc10 cluster_q7_pc10 0.2620445 7 10
# cluster_q9_pc10 cluster_q9_pc10 0.2152346 9 10


ggplot(brain_subset, aes(x = x, y = y, color = as.factor(cluster_q7_pc3))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=7")

#-------------------------------------------------------------------------------
##### k-means#####
set.seed(42)
load("Brain/DFPLC_151510_data_for_iIMPACT.RData")

clusters <- c(5, 7, 9)

# Preallocate lists to store results
results_q3 <- list()
results_q10 <- list()

for (c in clusters) {
  kmeans_res_pc3 <- kmeans(Y[, 1:3], c)

  # Store clustering labels
  df_q3 <- data.frame(
    spatial_cluster = kmeans_res_pc3$cluster,
    PCs = 3,
    q = c
  )
  results_q3[[as.character(c)]] <- df_q3

  kmeans_res_pc10 <- kmeans(Y[, 1:10], c)

  df_q10 <- data.frame(
    spatial_cluster = kmeans_res_pc10$cluster,
    PCs = 10,
    q = c
  )
  results_q10[[as.character(c)]] <- df_q10
}

z = ground_truth_df$type

brain_kmeans_df <- data.frame(loc,
                              pc3_k5 = results_q3[["5"]][["spatial_cluster"]],
                              pc3_k7 = results_q3[["7"]][["spatial_cluster"]],
                              pc3_k9 = results_q3[["9"]][["spatial_cluster"]],
                              pc10_k5 = results_q10[["5"]][["spatial_cluster"]],
                              pc10_k7 = results_q10[["7"]][["spatial_cluster"]],
                              pc10_k9 = results_q10[["9"]][["spatial_cluster"]],
                              z)

saveRDS(brain_kmeans_df, file = "Brain/kmeans_Brain_clustering_results.rds")

brain_subset <- subset(brain_kmeans_df, !is.na(z))
head(brain_subset)

clustering_cols <- grep("^pc", names(brain_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  round(adjustedRandIndex(brain_subset$z, brain_subset[[col]]), 3)
})

# Make dataframe and extract PC + K info
ari_df <- data.frame(
  Clustering = names(ari_scores),
  ARI = ari_scores
) %>%
  mutate(
    PC = str_extract(Clustering, "pc\\d+"),
    K = as.integer(str_extract(Clustering, "(?<=k)\\d+"))
  ) %>%
  arrange(PC, K)

# Split into two tables
ari_pc3 <- ari_df %>% filter(PC == "pc3") %>% arrange(desc(ARI))
ari_pc10 <- ari_df %>% filter(PC == "pc10") %>% arrange(desc(ARI))

# View results
print(ari_pc3)
# Clustering   ARI  PC K
# pc3_k5     pc3_k5 0.264 pc3 5
# pc3_k7     pc3_k7 0.253 pc3 7
# pc3_k9     pc3_k9 0.218 pc3 9

print(ari_pc10)
# Clustering   ARI   PC K
# pc10_k5    pc10_k5 0.257 pc10 5
# pc10_k9    pc10_k9 0.221 pc10 9
# pc10_k7    pc10_k7 0.211 pc10 7

ggplot(brain_subset, aes(x = x, y = y, color = as.factor(pc3_k7))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=7")


