import sys
import scanpy as sc
import numpy as np
import pandas as pd
import torch
import copy
import matplotlib.pyplot as plt
import warnings
import os
os.environ['R_LIBS_USER'] = "/scratch/user/varogovchenko/Rlibs"
sys.path.append("../STAMarker")
from stamarker.dataset import SpatialDataModule
from stamarker.pipeline import STAMarker, make_spatial_data
from stamarker.utils import parse_args
warnings.filterwarnings("ignore")

# --- 1. Set up logging ---
logfile = open("STAMarker_Lung_Xenium_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

# --- 2. Load and preprocess data ---
H5AD = "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/Data_VUILD96MF/Lung_xenium_processed.h5ad"
ann_data = sc.read_h5ad(H5AD)
ann_data.var_names_make_unique()
print(ann_data)

# Set up spatial coordinates
ann_data.obsm['spatial'] = ann_data.obs[['x_centroid', 'y_centroid']].values

# Use pre-computed highly variable genes
ann_data = ann_data[:, ann_data.var.highly_variable].copy()
data_module = make_spatial_data(ann_data)

# Data is already processed, but we still need to build the spatial graph
data_module.prepare_data(rad_cutoff=15, n_top_genes=None, min_cells=0) # n_top_genes=None to skip HVG selection

# --- 3. Configure and initialize the model ---
config = dict()
config.update(parse_args("../STAMarker/tutorial/_params/model.yaml"))
config.update(parse_args("../STAMarker/tutorial/_params/trainer.yaml"))
config["stagate"]["params"]["in_features"] = data_module.ann_data.shape[1]

# Disable GPU if not available
if not torch.cuda.is_available():
    config["stagate_trainer"]["gpus"] = None
    config["classifier_trainer"]["gpus"] = None

# Initialize the STAMarker model
model = STAMarker(5, "Lung_Xenium_STAMarker_output/", config)

# --- 4. Run the STAMarker pipeline ---
# Train autoencoders
model.train_auto_encoders(data_module)

# Perform clustering
n_clusters = 6
model.clustering(data_module, "mclust", n_clusters)
print("Finished m-clustering, starting consensus clustering...")
consensus_labels = model.consensus_clustering(n_clusters, show_plot=False, large_dataset=True)

# Train classifiers
model.train_classifiers(data_module, n_clusters, consensus_labels_path="consensus_labels.npy")

# Compute spatial maps
smaps = model.compute_smaps(data_module, return_recon=False)

# --- 5. Prepare data for validation ---
consensus_labels = np.load(model.save_dir + "/consensus_labels.npy")
ann_data = copy.copy(data_module.ann_data)
ann_data.obs["Consensus clustering"] = consensus_labels.astype(str)
n_class = np.max(consensus_labels) + 1
print("Num of spatial domains", n_class)

# --- 6. ARI calculation and comparison plotting ---
from sklearn.metrics import adjusted_rand_score

# Get ground truth labels from the AnnData object
true_labels = ann_data.obs["Annotation_Type"].dropna()
consensus_for_ari = ann_data.obs["Consensus clustering"][true_labels.index]

ari = adjusted_rand_score(true_labels, consensus_for_ari)
print(f"Adjusted Rand Index (ARI): {ari}")

# Generate and save comparison plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.embedding(ann_data, basis="spatial", s=25, color="Consensus clustering", ax=axes[0], show=False, title="Consensus Clustering")

# Create a temporary AnnData object for plotting true labels
temp_adata = ann_data[true_labels.index].copy()
sc.pl.embedding(temp_adata, basis="spatial", s=25, color="Annotation_Type", ax=axes[1], show=False, title="True Labels")

plt.tight_layout()
plt.savefig("Lung_Xenium_STAMarker.png")


# --- 7. Generate final output CSV ---
# Create a dataframe with coordinates, consensus labels, and true labels
output_df = pd.DataFrame(
    ann_data.obsm['spatial'],
    columns=['x', 'y'],
    index=ann_data.obs.index
)
output_df['barcode'] = output_df.index
output_df['STAMarker_label'] = ann_data.obs["Consensus clustering"]
if "Annotation_Type" in ann_data.obs.columns:
    output_df['true_label'] = ann_data.obs["Annotation_Type"]
else:
    output_df['true_label'] = np.nan

output_df.to_csv("STAMarker_Lung_Xenium_joint.csv", index=False)