# Load necessary libraries for single-cell and spatial transcriptomics analysis
library(scCustomize)          
library(SingleCellExperiment) 
library(data.table)           
library(ggplot2)              
library(scater)               
library(dplyr)

set.seed(101) # Set a random seed for reproducibility

# Define the directory where the raw data is stored
data_dir <- './Gut/Space_Ranger_Data_Gut'

# Read in the expression matrices from 10X hdf5 format
expression_matrices <- Read10X_h5_GEO(data_dir = data_dir)

# Create a Seurat object from the first matrix in the list
seurat_object = CreateSeuratObject(counts = expression_matrices[[1]])
seurat_object

# Convert the Seurat object to a SingleCellExperiment object
sce <- as.SingleCellExperiment(seurat_object)

# Log-transform the raw counts
log_count <- log(sce@assays@data@listData$counts + 1) 
# Convert the log-transformed matrix to a sparse matrix format (saves memory)
M1 <- as(log_count, "CsparseMatrix") #"dgCMatrix"

# Store the log-transformed counts back in the SingleCellExperiment object
sce@assays@data@listData$logcounts <- M1
rm(log_count) # Remove log_count to free up memory

# Use scran to model gene variability and identify highly variable genes (HVGs)
dec <- scran::modelGeneVar(sce)
# Get the top 2000 most variable genes for dimensionality reduction
top <- scran::getTopHVGs(dec, n = 2000)


# #-----Get the top 100 highly variable genes for visualization in Figure 1-----#
# top_plot <- scran::getTopHVGs(dec, n = 100)
# # Plot a heatmap of the top 100 highly variable genes
# scater::plotHeatmap(sce, 
#                     features = top_plot, 
#                     exprs_values = "logcounts", 
#                     color = viridis::viridis(100))
# #-----


set.seed(102)
# Perform PCA using the top 2000 highly variable genes, reducing to 15 components
sce <- scater::runPCA(sce, subset_row=top, ncomponents = 15)

# Extract loadings
# loadings <- attr(reducedDim(sce, "PCA"), "rotation")
# loadings_df <- as.data.frame(loadings)
# write.csv(loadings_df, "./Gut/PCA_loadings.csv", row.names = TRUE)


# Add BayesSpace metadata
sce <- BayesSpace::spatialPreprocess(sce, platform="Visium", skip.PCA=TRUE)

# Extract the PCA matrix from the SingleCellExperiment object
pca_reduced <- sce@int_colData@listData$reducedDims@listData
PCA <- pca_reduced[["PCA"]]

# Convert PCA matrix to a data frame
PCA_df <- data.frame(PCA)
# Add row names as a new column
PCA_df$ID <- rownames(PCA_df)

# Load the tissue position data from a CSV file
dt = fread("./Real Data/GSM5213484_V19S23-097_B1_S2_tissue_positions_list.csv.gz")

# Convert dt to a data frame
dt_df <- data.frame(dt)
colnames(dt_df)[1] <- "ID"

# Merge the PCA data frame with the tissue position data using the "ID" column
merged_data <- merge(PCA_df, dt_df, by = "ID")
head(merged_data)
rownames(merged_data) <- merged_data$ID # Set the row names of the merged data to the "ID" column

# Select only relevant columns (PCA components, x and y coordinates) from merged data
merged_data <- merged_data[, c(2:4, 8, 9)]
colnames(merged_data)[c(4,5)] <- c("x", "y") # Rename the x and y coordinate columns for clarity
head(merged_data)

### Plot heatmaps of the first three PCs
PC1 <- ggplot() + 
  geom_point(aes(x = x, y = y, col = PC1), data = merged_data) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") + 
  ggtitle('First PCA') +
  theme(legend.position = "none")
PC1

PC2 <- ggplot() + 
  geom_point(aes(x = x, y = y, col = PC2), data = merged_data) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") + 
  ggtitle('Second PCA') +
  theme(legend.position = "none")
PC2

PC3 <- ggplot() + 
  geom_point(aes(x = x, y = y, col = PC3), data = merged_data) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") + 
  ggtitle('Third PCA') +
  theme(legend.position = "none")
PC3


##### Manually remove muscle layer and add true hystological labels #####

# Manual selection of muscle cells for deletion was performed using the identify() 
# function. By plotting the spatial coordinates (x, y) of the dataset and 
# interactively selecting points corresponding to muscle cells, we obtained the 
# indices (out_delete_index). These indices were then used to remove the muscle 
# cells from the dataset, ensuring that subsequent analyses focused on relevant 
# tissue compartments.

plot(merged_data$x, merged_data$y)

out_delete_index <- identify(merged_data$x, merged_data$y)
out_delete_index

swiss_roll_wt_muscle_finaltouches1 = merged_data[-out_delete_index,]

# Then, we manually annotated all cells according to the provided histological 
# annotations, storing the results in the z_man variable for further analysis.

# Group the manual annotations into z variable
swiss_roll_wt_muscle_finaltouches1 <- swiss_roll_wt_muscle_finaltouches1 %>%
  mutate(z = case_when(
    z_man %in% c("Normal", "Normal with occasional alterations", "Damaged") ~ "Normal phenotype",
    z_man %in% c("Edema and inflammation", "Squamous area, mild edema in deeper layers") ~ "Edema phenotype",
    z_man %in% c("Inflammation and hyperplasia", "Occasional inflammation and hyperplasia", "Mild Inflammation") ~ "Inflammation phenotype",
    z_man == "Large lymphoid patch" ~ "Lymphoid phenotype",
    z_man == "Mild hyperchromasy and regional loss of goblet cells" ~ "Loss of goblet cells",
    TRUE ~ NA_character_  # This will handle any categories not explicitly listed
  ))

## Save the final dataset
# save(swiss_roll_wt_muscle_finaltouches1, 
#      file = "./Real Data//swiss_roll_wt_muscle_finaltouches1.RData")

load("./Real Data/swiss_roll_wt_muscle_finaltouches1.RData")
head(swiss_roll_wt_muscle_finaltouches1)

# The final dataset consists of the first three principal components (PC1, PC2, PC3), 
# which represent the reduced dimensionality of the gene expression data following PCA. 
# Each cell is assigned spatial coordinates ('x', 'y'), indicating its location within 
# the tissue section. The 'z_man' column contains manually assigned phenotypic labels 
# based on histological annotations provided in the original paper. 
# The 'z' column provides further grouped descriptions of these annotations used
# in our analysis.
# Overall, there are 2,568 spots left for the analysis.


# Scatter plot of cells by spatial coordinates for original hystological annotation
ggplot(swiss_roll_wt_muscle_finaltouches1, aes(x = x, y = y, color = z_man)) +
  geom_point() +  # Add points
  labs(title = "Scatter Plot with Color Coding by z_man",
       x = "x Coordinate",
       y = "y Coordinate") +
  theme_void() 

# Scatter plot of cells by spatial coordinates for grouped hystological annotation
ggplot(swiss_roll_wt_muscle_finaltouches1, aes(x = x, y = y, color = z)) +
  geom_point() +  # Add points
  labs(title = "Scatter Plot with Color Coding by z_man",
       x = "x Coordinate",
       y = "y Coordinate") +
  theme_void() 


##### Manually create the boundary #####

# Plot spatial coordinates and manually draw a boundary using the locator tool.
# The boundary is saved as a list of x and y coordinates and stored in an RData file.
plot(swiss_roll_wt_muscle_finaltouches1$x, swiss_roll_wt_muscle_finaltouches1$y)

boundary <- locator(n = 5000, type = "l")

boundary = list()
boundary$x = c(boundary$x)
boundary$y = c(boundary$y)

# save(boundary, file = "./Real Data/swiss_roll_wt_muscle_boundary.RData")
load("./Real Data/swiss_roll_wt_muscle_boundary.RData")
