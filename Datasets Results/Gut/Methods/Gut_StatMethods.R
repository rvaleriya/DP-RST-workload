##### COMPARE WITH OTHER METHODS #####
library(ggplot2)
library(dplyr)
library(purrr)
library(mclust)
library(tidyr)
library(stringr)
library(scCustomize)
#-------------------------------------------------------------------------------
set.seed(42)
### BayesSpace ###
library(BayesSpace)

data_dir <- 'Gut'
expression_matrices <- Read10X_h5_GEO(data_dir = data_dir)

load("Gut/swiss_roll_wt_muscle_finaltouches1.RData")
loc = swiss_roll_wt_muscle_finaltouches1[, 4:5]
names(loc) <- c("row", "col")

# Step 1: Identify spots (cells) with valid coordinates
valid_spots <- rownames(loc)[!is.na(loc$row) & !is.na(loc$col)]

# Step 2: Subset the expression matrix
expression_matrix_filtered <- expression_matrices[[1]][, valid_spots]

# Step 3: Check dimensions to confirm
dim(expression_matrix_filtered)

# To create object from single file
seurat_object = CreateSeuratObject(counts = expression_matrix_filtered, meta.data = loc)

# Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(seurat_object)

log_count <- log(sce@assays@data@listData$counts + 1)
M1 <- as(log_count, "CsparseMatrix") #"dgCMatrix"

sce@assays@data@listData$logcounts <- M1
rm(log_count)

dec <- scran::modelGeneVar(sce)
top <- scran::getTopHVGs(dec, n = 2000)

sce <- scater::runPCA(sce, subset_row=top, ncomponents = 15)

## Add BayesSpace metadata
sce <- BayesSpace::spatialPreprocess(sce, platform="Visium", skip.PCA=TRUE)

pca_reduced <- sce@int_colData@listData$reducedDims@listData


clusters <- c(3, 5, 7)

# Preallocate lists to store results
results_q3 <- list()
results_q10 <- list()

for (c in clusters) {
  # Clustering with q=3
  gut_q3 <- spatialCluster(sce, q=c, d=3, platform='Visium',
                              init.method="mclust", model="t", gamma=3,
                              nrep=50000, burn.in=40000, save.chain=FALSE)

  # Store clustering labels
  df_q3 <- data.frame(
    imagerow = gut_q3$row,
    imagecol = gut_q3$col,
    spatial_cluster = gut_q3$spatial.cluster,
    PCs = 3,
    q = c
  )
  results_q3[[as.character(c)]] <- df_q3

  # Clustering with q=10
  gut_q10 <- spatialCluster(sce, q=c, d=10, platform='Visium',
                               init.method="mclust", model="t", gamma=3,
                               nrep=50000, burn.in=40000, save.chain=FALSE)

  # Store clustering labels
  df_q10 <- data.frame(
    imagerow = gut_q10$row,
    imagecol = gut_q10$col,
    spatial_cluster = gut_q10$spatial.cluster,
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
saveRDS(merged_df, file = "Gut/BayesSpace_Gut_clustering_results.rds")


merged_df_renamed <- merged_df %>%
  rename(x = imagerow, y = imagecol)

loc = swiss_roll_wt_muscle_finaltouches1[, 4:5]
ground_truth_df <- cbind(loc, z = swiss_roll_wt_muscle_finaltouches1$z)

gut_merged_df <- merge(ground_truth_df, merged_df_renamed,
                         by = c("x", "y"),
                         all = TRUE)


gut_subset <- subset(gut_merged_df, !is.na(gut_merged_df$z))
head(gut_subset)

clustering_cols <- grep("^cluster_", names(gut_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  adjustedRandIndex(gut_subset$z, gut_subset[[col]])
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
# Clustering        ARI q PC
# cluster_q5_pc3 cluster_q5_pc3 0.12803357 5  3
# cluster_q7_pc3 cluster_q7_pc3 0.11164551 7  3
# cluster_q3_pc3 cluster_q3_pc3 0.01439407 3  3

print(ari_pc10)
# Clustering        ARI q PC
# cluster_q7_pc10 cluster_q7_pc10 0.08373477 7 10
# cluster_q5_pc10 cluster_q5_pc10 0.06458976 5 10
# cluster_q3_pc10 cluster_q3_pc10 0.04077333 3 10


ggplot(gut_subset, aes(x = x, y = y, color = as.factor(cluster_q5_pc3))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=7")


#-------------------------------------------------------------------------------
### SC.MEB ###
set.seed(42)

library(SC.MEB)
gut_df_wt_muscle <- readRDS("~/Desktop/DP.RST/Gut/gut_df_wt_muscle.rds")

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
Adj_sp <- getneighborhood_fast(as.matrix(gut_df_wt_muscle[, c("x", "y")]), cutoff = 1.2)

#### Run the SC-MEB in parallel
gut_SC.MEB_3pc = SC.MEB(as.matrix(gut_df_wt_muscle[, 1:3]), Adj_sp, beta_grid = beta_grid, K_set = K_set,
                             parallel = parallel, num_core = num_core, PX = PX,
                             maxIter_ICM = maxIter_ICM, maxIter = maxIter)
str(gut_SC.MEB_3pc[,1])
str(gut_SC.MEB_3pc[,2])
str(gut_SC.MEB_3pc[,3])

gut_SC.MEB_10pc = SC.MEB(as.matrix(gut_df_wt_muscle[, 1:10]), Adj_sp, beta_grid = beta_grid, K_set = K_set,
                              parallel = parallel, num_core = num_core, PX = PX,
                              maxIter_ICM = maxIter_ICM, maxIter = maxIter)

# The item 'x' is clustering label.
# The item 'ell' is the opposite log-likelihood for each beta and K.
# The item 'mu' is the mean of each component.
# The item 'sigma' is the variance of each component.
# The item 'gam' is the posterior probability.
# The item 'beta' is the estimated smoothing parameter.


# #### Selecting the number of clusters using BIC
# selectKPlot(gut_fit, K_set = K_set, criterion = "BIC")
#
# #### Selecting the number of clusters using Modified BIC
# selectKPlot(gut_fit, K_set = K_set, criterion = "MBIC")
#
# out = selectK(gut_fit, K_set = K_set, criterion = "BIC")
# ClusterPlot(out, loc)
#
# gut_SC.MEB_df <- data.frame(loc, out[["best_K_label"]], z)
gut_SC.MEB_df <- data.frame(as.matrix(gut_df_wt_muscle[, c("x", "y")]),
                                 pc3_k3 = gut_SC.MEB_3pc[,1]$x,
                                 pc3_k5 = gut_SC.MEB_3pc[,2]$x,
                                 pc3_k7 = gut_SC.MEB_3pc[,3]$x,
                                 pc10_k3 = gut_SC.MEB_10pc[,1]$x,
                                 pc10_k5 = gut_SC.MEB_10pc[,2]$x,
                                 pc10_k7 = gut_SC.MEB_10pc[,3]$x,
                            z = gut_df_wt_muscle$z)
head(gut_SC.MEB_df)
saveRDS(gut_SC.MEB_df, file = "Gut/SC-MEB_Gut_clustering_results.rds")

gut_subset <- subset(gut_SC.MEB_df, !is.na(z))
head(gut_subset)

clustering_cols <- grep("^pc", names(gut_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  round(adjustedRandIndex(gut_subset$z, gut_subset[[col]]), 3)
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
# pc3_k3     pc3_k3 0.169 pc3 3
# pc3_k5     pc3_k5 0.112 pc3 5
# pc3_k7     pc3_k7 0.100 pc3 7

print(ari_pc10)
# Clustering   ARI   PC K
# pc10_k5    pc10_k5 0.110 pc10 5
# pc10_k7    pc10_k7 0.109 pc10 7
# pc10_k3    pc10_k3 0.094 pc10 3


ggplot(gut_subset, aes(x = x, y = y, color = as.factor(pc3_k5))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=5")

#-------------------------------------------------------------------------------
### DR.SC ###
set.seed(42)

library(Seurat)
library(DR.SC)
### Create Seurat object ###
# Create count matrix
data_dir <- 'Gut/Space_Ranger_Data_Gut'
expression_matrix <- Read10X_h5_GEO(data_dir = data_dir)

load("Gut/swiss_roll_wt_muscle_finaltouches1.RData")
loc = swiss_roll_wt_muscle_finaltouches1[, 4:5]

colnames(loc) <- c("row", "col")

# Create Seurat object
seurat_object = CreateSeuratObject(counts = expression_matrix, meta.data = loc)
head(seurat_object)

# First, extract the metadata
meta <- seurat_object@meta.data
# Filter metadata to exclude rows with NA in 'row' or 'col'
filtered_meta <- meta[!is.na(meta$row) & !is.na(meta$col), ]
# Subset Seurat object to keep only cells (spots) with valid coordinates
seurat_object_filtered <- subset(seurat_object, cells = rownames(filtered_meta))
head(seurat_object_filtered)
dim(seurat_object_filtered)

# standard log-normalization
Swiss <- NormalizeData(seurat_object_filtered, verbose = F)
# choose 2000 highly variable features
seu <- FindVariableFeatures(Swiss, nfeatures = 2000, verbose = F)


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
saveRDS(merged_df, file = "Gut/DR-SC_Gut_clustering_results.rds")

merged_df_renamed <- merged_df %>%
  rename(x = imagerow, y = imagecol)

ground_truth_df <- swiss_roll_wt_muscle_finaltouches1[, c("x", "y", "z")]

gut_merged_df <- merge(ground_truth_df, merged_df_renamed,
                          by = c("x", "y"),
                          all = TRUE)

gut_subset <- subset(gut_merged_df, !is.na(gut_merged_df$z))
head(gut_subset)

clustering_cols <- grep("^cluster_", names(gut_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  adjustedRandIndex(gut_subset$z, gut_subset[[col]])
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
# cluster_q3_pc3 cluster_q3_pc3 0.1648175 3  3
# cluster_q5_pc3 cluster_q5_pc3 0.1530579 5  3
# cluster_q7_pc3 cluster_q7_pc3 0.1449124 7  3

print(ari_pc10)
# Clustering       ARI q PC
# cluster_q3_pc10 cluster_q3_pc10 0.2713140 3 10
# cluster_q5_pc10 cluster_q5_pc10 0.1950344 5 10
# cluster_q7_pc10 cluster_q7_pc10 0.1320986 7 10

ggplot(gut_subset, aes(x = x, y = y, color = as.factor(cluster_q5_pc3))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=3")


#-------------------------------------------------------------------------------
##### k-means#####
set.seed(42)
gut_df_wt_muscle <- readRDS("~/Desktop/DP.RST/Gut/gut_df_wt_muscle.rds")

loc = gut_df_wt_muscle[, c("x", "y")]

clusters <- c(3, 5, 7)

# Preallocate lists to store results
results_q3 <- list()
results_q10 <- list()

for (c in clusters) {
  kmeans_res_pc3 <- kmeans(gut_df_wt_muscle[, 1:3], c)

  # Store clustering labels
  df_q3 <- data.frame(
    spatial_cluster = kmeans_res_pc3$cluster,
    PCs = 3,
    q = c
  )
  results_q3[[as.character(c)]] <- df_q3

  kmeans_res_pc10 <- kmeans(gut_df_wt_muscle[, 1:10], c)

  df_q10 <- data.frame(
    spatial_cluster = kmeans_res_pc10$cluster,
    PCs = 10,
    q = c
  )
  results_q10[[as.character(c)]] <- df_q10
}

gut_kmeans_df <- data.frame(loc,
                               pc3_k3 = results_q3[["3"]][["spatial_cluster"]],
                               pc3_k5 = results_q3[["5"]][["spatial_cluster"]],
                               pc3_k7 = results_q3[["7"]][["spatial_cluster"]],
                               pc10_k3 = results_q10[["3"]][["spatial_cluster"]],
                               pc10_k5 = results_q10[["5"]][["spatial_cluster"]],
                               pc10_k7 = results_q10[["7"]][["spatial_cluster"]],
                               z = gut_df_wt_muscle$z)

saveRDS(gut_kmeans_df, file = "Gut/kmeans_Gut_clustering_results.rds")

gut_subset <- subset(gut_kmeans_df, !is.na(z))
head(gut_subset)

clustering_cols <- grep("^pc", names(gut_subset), value = TRUE)

ari_scores <- sapply(clustering_cols, function(col) {
  round(adjustedRandIndex(gut_subset$z, gut_subset[[col]]), 3)
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
# Clustering    ARI  PC K
# pc3_k5     pc3_k5  0.030 pc3 5
# pc3_k7     pc3_k7  0.030 pc3 7
# pc3_k3     pc3_k3 -0.001 pc3 3

print(ari_pc10)
# Clustering   ARI   PC K
# pc10_k5    pc10_k5 0.023 pc10 5
# pc10_k7    pc10_k7 0.022 pc10 7
# pc10_k3    pc10_k3 0.004 pc10 3

ggplot(gut_subset, aes(x = x, y = y, color = as.factor(pc3_k5))) +
  geom_point(size = 1) +
  scale_color_discrete(name = "Cluster") +
  coord_fixed() +  # keep aspect ratio
  theme_minimal() +
  ggtitle("Clustering result: PC=3, K=7")








