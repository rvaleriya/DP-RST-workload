# Load necessary libraries
library(SingleCellExperiment)
library("DR.SC")
library(Seurat)
library(ggplot2)
library(mclust)
library(paletteer)
library(dplyr)
library(rprojroot)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

# Set seed
set.seed(9362)

# Load custom functions
source("./Codes/Custom Functions/ComplexDomainFun.R")
source("./Codes/Custom Functions/plotting_fun.R")

### Load Data ###
load("./Simulations/Swiss_Roll_3k/Swiss_Roll_sim_3k_data_1.RData")
load("./Simulations/Swiss_Roll_3k/SR_sim_PT_30reps.RData")

# Get seeds from the BAST run for initialization
seeds <- c()
for (i in 1:length(SR_sim_PT_reps)){
  seeds[i] <- SR_sim_PT_reps[[i]][["seed"]]
}
rm(SR_sim_PT_reps) # Remove BASTinit data

# Extract first three PCAs
Y_sample = sim_data[, 5:7]
loc = sim_data[, 1:2]

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
SR_sim_DR.SC_reps <- list()

for (s in 1:length(seeds)){
  
  set.seed(seeds[s])
  
  DR.SC_fit <- DR.SC_fit(Y_std, K=3, Adj_sp=Adj_sp, q=3,
                         error.heter= TRUE, beta_grid=seq(0.5, 5, by=0.5),
                         maxIter=25, epsLogLik=1e-5, verbose=FALSE, maxIter_ICM=6,
                         wpca.int=FALSE, int.model="EEE", approxPCA=FALSE, coreNum = 5)
  
  SR_sim_DR.SC_reps[[s]] <- DR.SC_fit
}

# Save the results of the repetitive runs
save(SR_sim_DR.SC_reps, file = "./Simulations/Swiss_Roll_3k/SR_sim_DR.SC_30reps.RData")
