##### Human prostate Cancer: Ductal Carcinoma In Situ, Invasive Carcinoma (FFPE) #####
# https://www.10xgenomics.com/datasets/human-prostate-cancer-ductal-carcinoma-in-situ-invasive-carcinoma-ffpe-1-standard-1-3-0
# https://www.10xgenomics.com/datasets/human-prostate-cancer-adenocarcinoma-with-invasive-carcinoma-ffpe-1-standard-1-3-0

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
sce_prostate <- readVisium("Prostate/Space_Ranger_Data_Prostate")

prostate <- spatialPreprocess(sce_prostate)

# prostate_tune <- qTune(prostate, qs=seq(2, 20), d=3)
# qPlot(prostate_tune)

clusters <- c(3, 5, 7)

# Preallocate lists to store results
results_q3 <- list()
results_q10 <- list()

for (c in clusters) {
  # Clustering with q=3
  prostate_q3 <- spatialCluster(prostate, q=c, d=3, platform='Visium',
                             init.method="mclust", model="t", gamma=3,
                             nrep=50000, burn.in=40000, save.chain=FALSE)

  # Store clustering labels
  df_q3 <- data.frame(
    imagerow = prostate_q3$imagerow,
    imagecol = prostate_q3$imagecol,
    spatial_cluster = prostate_q3$spatial.cluster,
    PCs = 3,
    q = c
  )
  results_q3[[as.character(c)]] <- df_q3

  # Clustering with q=10
  prostate_q10 <- spatialCluster(prostate, q=c, d=10, platform='Visium',
                              init.method="mclust", model="t", gamma=3,
                              nrep=50000, burn.in=40000, save.chain=FALSE)

  # Store clustering labels
  df_q10 <- data.frame(
    imagerow = prostate_q10$imagerow,
    imagecol = prostate_q10$imagecol,
    spatial_cluster = prostate_q10$spatial.cluster,
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
saveRDS(merged_df, file = "Prostate/BayesSpace_Prostate_clustering_results.rds")


load("Prostate/10x_prostate_cancer_iIMPACT_data.RData")
prostate_true_labels <- data.frame(loc, z)


merged_df_renamed <- merged_df %>%
  rename(x = imagerow, y = imagecol)

prostate_merged_df <- merge(prostate_true_labels, merged_df_renamed,
                         by = c("x", "y"),
                         all = TRUE)

prostate_subset <- subset(prostate_merged_df, !is.na(prostate_merged_df$z))
head(prostate_subset)

clustering_cols <- grep("^cluster_", names(prostate_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  adjustedRandIndex(prostate_subset$z, prostate_subset[[col]])
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
# cluster_q5_pc3 cluster_q5_pc3 0.6759571 5  3
# cluster_q7_pc3 cluster_q7_pc3 0.6055482 7  3
# cluster_q3_pc3 cluster_q3_pc3 0.5953267 3  3

print(ari_pc10)

# Clustering       ARI q PC
# cluster_q5_pc10 cluster_q5_pc10 0.6277993 5 10
# cluster_q7_pc10 cluster_q7_pc10 0.5837904 7 10
# cluster_q3_pc10 cluster_q3_pc10 0.5405020 3 10

ggplot(prostate_subset, aes(x = x, y = y, color = as.factor(cluster_q3_pc3))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=3")



##### SC-MEB #####
set.seed(42)

library(SC.MEB)
load("Prostate/10x_prostate_cancer_iIMPACT_data.RData")

### Set the parameters
#platform = "ST"
beta_grid = seq(0, 5, 0.2) #Vector specifying the smoothness of Random Markov Field
K_set= c(3, 5, 7) # 2:10 #Numbers of mixture components
parallel=TRUE #Logical value specifing the run the model in parallel
num_core = 3
PX = TRUE #Logical value for paramter expansion in EM algorithm
maxIter_ICM = 10 #Maximum iteration of ICM algorithm
maxIter = 50 #Maximum iteration of EM algorithm

#### Calculating the neighborhood from coordinates
Adj_sp <- getneighborhood_fast(as.matrix(loc), cutoff = 1.2)

#### Run the SC-MEB in parallel
prostate_SC.MEB_3pc = SC.MEB(Y[, 1:3], Adj_sp, beta_grid = beta_grid, K_set = K_set,
                          parallel = parallel, num_core = num_core, PX = PX,
                          maxIter_ICM = maxIter_ICM, maxIter = maxIter)
str(prostate_SC.MEB_3pc[,1])
str(prostate_SC.MEB_3pc[,2])
str(prostate_SC.MEB_3pc[,3])

prostate_SC.MEB_10pc = SC.MEB(Y[, 1:10], Adj_sp, beta_grid = beta_grid, K_set = K_set,
                           parallel = parallel, num_core = num_core, PX = PX,
                           maxIter_ICM = maxIter_ICM, maxIter = maxIter)

# The item 'x' is clustering label.
# The item 'ell' is the opposite log-likelihood for each beta and K.
# The item 'mu' is the mean of each component.
# The item 'sigma' is the variance of each component.
# The item 'gam' is the posterior probability.
# The item 'beta' is the estimated smoothing parameter.


# #### Selecting the number of clusters using BIC
# selectKPlot(prostate_fit, K_set = K_set, criterion = "BIC")
#
# #### Selecting the number of clusters using Modified BIC
# selectKPlot(prostate_fit, K_set = K_set, criterion = "MBIC")
#
# out = selectK(prostate_fit, K_set = K_set, criterion = "BIC")
# ClusterPlot(out, loc)
#
# prostate_SC.MEB_df <- data.frame(loc, out[["best_K_label"]], z)
prostate_SC.MEB_df <- data.frame(loc,
                              pc3_k3 = prostate_SC.MEB_3pc[,1]$x,
                              pc3_k5 = prostate_SC.MEB_3pc[,2]$x,
                              pc3_k7 = prostate_SC.MEB_3pc[,3]$x,
                              pc10_k3 = prostate_SC.MEB_10pc[,1]$x,
                              pc10_k5 = prostate_SC.MEB_10pc[,2]$x,
                              pc10_k7 = prostate_SC.MEB_10pc[,3]$x,
                              z)
head(prostate_SC.MEB_df)
saveRDS(prostate_SC.MEB_df, file = "Prostate/SC-MEB_Prostate_clustering_results.rds")

prostate_subset <- subset(prostate_SC.MEB_df, !is.na(z))
head(prostate_subset)

clustering_cols <- grep("^pc", names(prostate_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  round(adjustedRandIndex(prostate_subset$z, prostate_subset[[col]]), 3)
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
# pc3_k3     pc3_k3 0.764 pc3 3
# pc3_k5     pc3_k5 0.540 pc3 5
# pc3_k7     pc3_k7 0.460 pc3 7

print(ari_pc10)
# Clustering   ARI   PC K
# pc10_k3    pc10_k3 0.754 pc10 3
# pc10_k5    pc10_k5 0.563 pc10 5
# pc10_k7    pc10_k7 0.508 pc10 7


ggplot(prostate_subset, aes(x = x, y = y, color = as.factor(pc3_k3))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=7")

##### DR-SC #####
set.seed(42)

library(Seurat)
library(DR.SC)
### Create Seurat object ###
# Create count matrix
data_dir_prostate_count <- 'Prostate/Space_Ranger_Data_prostate/filtered_feature_bc_matrix'
list.files(data_dir_prostate_count) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir_prostate_count)

# Load spatial data
tissue_positions_list <- readr::read_csv("Prostate/Space_Ranger_Data_prostate/spatial/tissue_positions_list.csv",
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
prostate <- NormalizeData(seurat_object, verbose = F)
# choose 2000 highly variable features
seu <- FindVariableFeatures(prostate, nfeatures = 2000, verbose = F)

# Fit DR-SC model using 2000 highly variable features
clusters <- c(3, 5, 7)

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
saveRDS(merged_df, file = "Prostate/DR-SC_Prostate_clustering_results.rds")

load("Prostate/10x_prostate_cancer_iIMPACT_data.RData")

merged_df_renamed <- merged_df %>%
  rename(x = imagerow, y = imagecol)

ground_truth_df <- cbind(loc, z)

prostate_merged_df <- merge(ground_truth_df, merged_df_renamed,
                         by = c("x", "y"),
                         all = TRUE)

prostate_subset <- subset(prostate_merged_df, !is.na(prostate_merged_df$z))
head(prostate_subset)

clustering_cols <- grep("^cluster_", names(prostate_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  adjustedRandIndex(prostate_subset$z, prostate_subset[[col]])
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
# cluster_q3_pc3 cluster_q3_pc3 0.7346611 3  3
# cluster_q5_pc3 cluster_q5_pc3 0.5644265 5  3
# cluster_q7_pc3 cluster_q7_pc3 0.3963085 7  3

print(ari_pc10)
# Clustering       ARI q PC
# cluster_q3_pc10 cluster_q3_pc10 0.6336794 3 10
# cluster_q5_pc10 cluster_q5_pc10 0.5307230 5 10
# cluster_q7_pc10 cluster_q7_pc10 0.5002655 7 10


ggplot(prostate_subset, aes(x = x, y = y, color = as.factor(cluster_q3_pc3))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=3")


#-------------------------------------------------------------------------------
##### k-means#####
set.seed(42)
load("Prostate/10x_prostate_cancer_iIMPACT_data.RData")

clusters <- c(3, 5, 7)

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

prostate_kmeans_df <- data.frame(loc,
                              pc3_k3 = results_q3[["3"]][["spatial_cluster"]],
                              pc3_k5 = results_q3[["5"]][["spatial_cluster"]],
                              pc3_k7 = results_q3[["7"]][["spatial_cluster"]],
                              pc10_k3 = results_q10[["3"]][["spatial_cluster"]],
                              pc10_k5 = results_q10[["5"]][["spatial_cluster"]],
                              pc10_k7 = results_q10[["7"]][["spatial_cluster"]],
                              z)

saveRDS(prostate_kmeans_df, file = "Prostate/kmeans_Prostate_clustering_results.rds")

prostate_subset <- subset(prostate_kmeans_df, !is.na(z))
head(prostate_subset)

clustering_cols <- grep("^pc", names(prostate_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  round(adjustedRandIndex(prostate_subset$z, prostate_subset[[col]]), 3)
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
# pc3_k3     pc3_k3 0.508 pc3 3
# pc3_k5     pc3_k5 0.379 pc3 5
# pc3_k7     pc3_k7 0.355 pc3 7

print(ari_pc10)
# Clustering   ARI   PC K
# pc10_k3    pc10_k3 0.509 pc10 3
# pc10_k5    pc10_k5 0.375 pc10 5
# pc10_k7    pc10_k7 0.354 pc10 7

ggplot(prostate_subset, aes(x = x, y = y, color = as.factor(pc3_k3))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=7")








