# Load necessary libraries

library(SingleCellExperiment)
library("DR.SC")
library(Seurat)
library(dplyr)

#-------------------------------------------------------------------------------

# Set working directory to the root of the R project
setwd("~/Desktop/DP-RST-workload")

# Set seed
set.seed(9362)

### Load Data ###
Curl_sim_data <- readr::read_csv("./New_Simulations/Curl/Curl_sim_data.csv")

#-------------------------------------------------------------------------------

seeds <- c(141, 549, 75, 492, 676, 179, 587, 592, 601, 916,
           518, 339, 921, 423, 330, 388, 273, 286, 61, 807,
           283, 127, 165, 952, 311, 597, 473, 594, 605, 101)

#-------------------------------------------------------------------------------
# Extract first three PCAs
Y_sample = Curl_sim_data[, 5:7] # 5:7 for PCs 1-3 and 5:14 for PCs 1-10
head(Y_sample)

loc = Curl_sim_data[, 1:2]
head(loc)

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

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

# Create dummy count data (e.g., one gene across all cells)
dummy_counts <- matrix(1, nrow = 2, ncol = ncol(sce))
rownames(dummy_counts) <- c("Dummy_Gene1", "Dummy_Gene2")
colnames(dummy_counts) <- colnames(sce)

colnames(coords) <- c("row", "col")

### Create Seurat object 
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

# Create a list to store results from different runs
Curl_sim_DR.SC_p3_reps <- list()

for (s in 1:length(seeds)){
  
  set.seed(seeds[s])
  
  DR.SC_fit <- DR.SC_fit(Y_std, K=10, Adj_sp=Adj_sp, q=p,
                         error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
                         maxIter=25, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,
                         wpca.int=FALSE, int.model="EEE", approxPCA=FALSE, coreNum = 5)
  
  Curl_sim_DR.SC_p3_reps[[s]] <- DR.SC_fit
}

# Save the results of the repetitive runs
save(Curl_sim_DR.SC_p3_reps, file = "./New_Simulations/Curl/Curl_results/Curl_sim_DR.SC_p3_30reps.RData")


#-------------------------------------------------------------------------------
# Extract first three PCAs
Y_sample = Curl_sim_data[, 5:14] # 5:7 for PCs 1-3 and 5:14 for PCs 1-10
head(Y_sample)

loc = Curl_sim_data[, 1:2]
head(loc)

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

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

# Create dummy count data (e.g., one gene across all cells)
dummy_counts <- matrix(1, nrow = 2, ncol = ncol(sce))
rownames(dummy_counts) <- c("Dummy_Gene1", "Dummy_Gene2")
colnames(dummy_counts) <- colnames(sce)

colnames(coords) <- c("row", "col")

### Create Seurat object 
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

# Create a list to store results from different runs
Curl_sim_DR.SC_p10_reps <- list()

for (s in 1:length(seeds)){
  
  set.seed(seeds[s])
  
  DR.SC_fit <- DR.SC_fit(Y_std, K=10, Adj_sp=Adj_sp, q=p,
                         error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
                         maxIter=25, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,
                         wpca.int=FALSE, int.model="EEE", approxPCA=FALSE, coreNum = 5)
  
  Curl_sim_DR.SC_p10_reps[[s]] <- DR.SC_fit
}

# Save the results of the repetitive runs
save(Curl_sim_DR.SC_p10_reps, file = "./New_Simulations/Curl/Curl_results/Curl_sim_DR.SC_p10_30reps.RData")

