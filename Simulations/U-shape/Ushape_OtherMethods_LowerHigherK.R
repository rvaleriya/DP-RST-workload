### Load Libraries ###
library(SingleCellExperiment)
library(BayesSpace)
library("DR.SC")
library("SC.MEB")
library(Seurat)
library(ggplot2)
library(dplyr)
library(mclust)
library(paletteer)

# Set seed
set.seed(9362)

### Load Data ###
load("/UshapeSim_coords_boundary1.RData")
load("/U_sim_PT_30reps.RData")

seeds <- c()
for (i in 1:1){
  seeds[i] <- U_sim_PT_reps[[i]][["seed"]]
}
rm(U_sim_PT_reps)

# Extract first three PCAs
Y_sample = sim_data[, 4:6]
loc = sim_data[, 1:2]

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

#-------------------------------------------------------------------------------
##### BAYES SPACE #####
rownames(Y_std) <- paste0("Cell", 1:n)
colnames(Y_std) <- paste0("PC", 1:p)
rownames(coords) <- rownames(Y_std)

# Create the SingleCellExperiment object
sce <- SingleCellExperiment(
  reducedDims = list(PCA = Y_std),  # Add PCA data as an assay
  colData = coords        # Add spatial coordinates as column data
)

# Rename the coordinates columns to match BayesSpace expectations
colnames(colData(sce)) <- c("col", "row")

set.seed(seeds[1])

BS_results_4c <- spatialCluster(sce, use.dimred = "PCA", d = 3, q = 4, platform = "ST")
BS_results_8c <- spatialCluster(sce, use.dimred = "PCA", d = 3, q = 8, platform = "ST")

BayesSpace_acc_4c <- adjustedRandIndex(sim_data$cluster, BS_results_4c$spatial.cluster)
BayesSpace_acc_4c # 0.3365667

BayesSpace_acc_8c <- adjustedRandIndex(sim_data$cluster, BS_results_8c$spatial.cluster)
BayesSpace_acc_8c # 0.276722

#-------------------------------------------------------------------------------
##### DR-SC #####
# Create dummy count data (e.g., one gene across all cells)
dummy_counts <- matrix(1, nrow = 2, ncol = ncol(sce))
rownames(dummy_counts) <- c("Dummy_Gene1", "Dummy_Gene2")
colnames(dummy_counts) <- colnames(sce)

colnames(coords) <- c("row", "col")

seurat_obj <- CreateSeuratObject(counts = dummy_counts)
# Add spatial coordinates to the metadata
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, coords)
# Verify that the coordinates are added
print(head(seurat_obj@meta.data))

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  select(row, col, everything())

# Verify the reordering
print(head(seurat_obj@meta.data))

# Create DimReduc object with PCA data
pca_object <- CreateDimReducObject(embeddings = Y_std, key = "PC_", 
                                   assay = DefaultAssay(seurat_obj))
# Add PCA data to Seurat object
seurat_obj[["pca"]] <- pca_object

Adj_sp <- getAdj(seurat_obj, platform = 'ST')

DR.SC_fit_4c <- DR.SC_fit(Y_std, K=4, Adj_sp=Adj_sp, q=3,
                       error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
                       maxIter=25, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,
                       wpca.int=FALSE, int.model="EEE", approxPCA=FALSE, coreNum = 5)


DR.SC_fit_8c <- DR.SC_fit(Y_std, K=8, Adj_sp=Adj_sp, q=3,
                          error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
                          maxIter=25, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,
                          wpca.int=FALSE, int.model="EEE", approxPCA=FALSE, coreNum = 5)


DR.SC_acc_4c <- adjustedRandIndex(sim_data$cluster, 
                                  DR.SC_fit_4c[["Objdrsc"]][[1]][["cluster"]])
DR.SC_acc_4c # 0.2274286

DR.SC_acc_8c <- adjustedRandIndex(sim_data$cluster, 
                                  DR.SC_fit_8c[["Objdrsc"]][[1]][["cluster"]])
DR.SC_acc_8c # 0.1732929

#-------------------------------------------------------------------------------
##### SC-MEB #####
### Set the parameters ###
beta_grid = seq(0, 5, 0.2) # Vector specifying the smoothness of Random Markov Field
parallel=TRUE #Logical value specifing the run the model in parallel
num_core = 3
PX = TRUE # Logical value for paramter expansion in EM algorithm
maxIter_ICM = 10 # Maximum iteration of ICM algorithm
maxIter = 50 # Maximum iteration of EM algorithm

#### Calculating the neighborhood from coordinates
Adj_sp  <- find_neighbors2(sce, platform = "ST")

SC.MEB_fit_4c = SC.MEB(y = Y_std, Adj_sp, beta_grid = beta_grid, K_set = 4, 
                    parallel = parallel, num_core = num_core, PX = PX, 
                    maxIter_ICM = maxIter_ICM, maxIter = maxIter)

SC.MEB_fit_8c = SC.MEB(y = Y_std, Adj_sp, beta_grid = beta_grid, K_set = 8, 
                       parallel = parallel, num_core = num_core, PX = PX, 
                       maxIter_ICM = maxIter_ICM, maxIter = maxIter)


SC.MEB_acc_4c <- adjustedRandIndex(sim_data$cluster, SC.MEB_fit_4c[[1]])
SC.MEB_acc_4c # 0.3192488

SC.MEB_acc_8c <- adjustedRandIndex(sim_data$cluster, SC.MEB_fit_8c[[1]])
SC.MEB_acc_8c # 0.255098

#-------------------------------------------------------------------------------
##### k-means #####

kmeans_fit_4c <- kmeans(Y_std, 4)$cluster

kmeans_fit_8c <- kmeans(Y_std, 8)$cluster

kmeans_acc_4c <- adjustedRandIndex(sim_data$cluster, kmeans_fit_4c)
kmeans_acc_4c # 0.3139693

kmeans_acc_8c <- adjustedRandIndex(sim_data$cluster, kmeans_fit_8c)
kmeans_acc_8c # 0.2414774






