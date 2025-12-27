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
from upsetplot import plot, from_contents
from scanpy.plotting.palettes import vega_20_scanpy
from stamarker.dataset import SpatialDataModule
from stamarker.pipeline import STAMarker, make_spatial_data
from stamarker.utils import parse_args, select_svgs
warnings.filterwarnings("ignore")

# --- 1. Set up logging ---
logfile = open("STAMarker_osmFISH_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

# --- 2. Load and preprocess data ---
H5AD = "/scratch/user/varogovchenko/BASTION_HPRC/osmFISH/osmfish_20251122011024.h5ad"
ann_data = sc.read_h5ad(H5AD)
ann_data.var_names_make_unique()
print(ann_data)

# Use adaptive number of genes based on dataset size
n_genes = ann_data.shape[1]
n_top_hvg = min(3000, n_genes)  # Don't exceed available genes
n_top_use = min(2000, n_genes)  # Don't exceed available genes

sc.pp.highly_variable_genes(ann_data, flavor="seurat_v3", n_top_genes=n_top_hvg)
# sc.pp.normalize_total(ann_data, target_sum=1e4)  # Already log normalized

# Set up spatial coordinates if not already present
if 'spatial' not in ann_data.obsm:
    if 'X' in ann_data.obs.columns and 'Y' in ann_data.obs.columns:
        ann_data.obsm['spatial'] = ann_data.obs[['X', 'Y']].values
    elif 'x' in ann_data.obs.columns and 'y' in ann_data.obs.columns:
        ann_data.obsm['spatial'] = ann_data.obs[['x', 'y']].values
    else:
        raise ValueError("Spatial coordinates (X/Y or x/y) not found in ann_data.obs")

data_module = make_spatial_data(ann_data)
data_module.prepare_data(rad_cutoff=150, n_top_genes=n_top_use, min_cells=0)

# --- 3. Configure and initialize the model ---
config = dict()
config.update(parse_args("../STAMarker/tutorial/_params/model.yaml"))
config.update(parse_args("../STAMarker/tutorial/_params/trainer.yaml"))
config["stagate"]["params"]["in_features"] = data_module.ann_data.shape[1]

# Set accelerator to cpu and remove conflicting GPU settings from the config file
config["stagate_trainer"]["accelerator"] = "cpu"
if "gpus" in config["stagate_trainer"]:
    del config["stagate_trainer"]["gpus"]
if "auto_select_gpus" in config["stagate_trainer"]:
    del config["stagate_trainer"]["auto_select_gpus"]

config["classifier_trainer"]["accelerator"] = "cpu"
if "gpus" in config["classifier_trainer"]:
    del config["classifier_trainer"]["gpus"]
if "auto_select_gpus" in config["classifier_trainer"]:
    del config["classifier_trainer"]["auto_select_gpus"]

# Initialize the STAMarker model
model = STAMarker(5, "osmFISH_STAMarker_output/", config)

# --- 4. Run the STAMarker pipeline ---
# Train autoencoders
model.train_auto_encoders(data_module)

# Perform clustering
n_clusters = 11
model.clustering(data_module, "mclust", n_clusters)
consensus_labels = model.consensus_clustering(n_clusters, show_plot=False)

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
if "ground_truth" in ann_data.obs.columns:
    from sklearn.metrics import adjusted_rand_score

    # Get ground truth labels from the AnnData object
    true_labels = ann_data.obs["ground_truth"].dropna()
    if not true_labels.empty:
        consensus_for_ari = ann_data.obs["Consensus clustering"][true_labels.index]

        ari = adjusted_rand_score(true_labels, consensus_for_ari)
        print(f"Adjusted Rand Index (ARI): {ari}")

        # Generate and save comparison plot
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.embedding(ann_data, basis="spatial", s=25, color="Consensus clustering", ax=axes[0], show=False, title="Consensus Clustering")

        # Create a temporary AnnData object for plotting true labels
        temp_adata = ann_data[true_labels.index].copy()
        sc.pl.embedding(temp_adata, basis="spatial", s=25, color="ground_truth", ax=axes[1], show=False, title="True Labels")

        plt.tight_layout()
        plt.savefig("osmFISH_STAMarker.png")
    else:
        plt.figure(figsize=(6, 5))
        sc.pl.embedding(ann_data, basis="spatial", s=25, color="Consensus clustering", show=False, title="Consensus Clustering")
        plt.savefig("osmFISH_STAMarker.png")
else:
    plt.figure(figsize=(6, 5))
    sc.pl.embedding(ann_data, basis="spatial", s=25, color="Consensus clustering", show=False, title="Consensus Clustering")
    plt.savefig("osmFISH_STAMarker.png")


# --- 7. Generate final output CSV ---
# Create a dataframe with coordinates, consensus labels, and true labels
output_df = pd.DataFrame(
    ann_data.obsm['spatial'],
    columns=['x', 'y'],
    index=ann_data.obs.index
)
output_df['barcode'] = output_df.index
output_df['STAMarker_label'] = ann_data.obs["Consensus clustering"]
if "ground_truth" in ann_data.obs.columns:
    output_df['true_label'] = ann_data.obs["ground_truth"]
else:
    output_df['true_label'] = np.nan

output_df.to_csv("STAMarker_osmFISH_joint.csv", index=False)
