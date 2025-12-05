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
# logfile = open("STAMarker_Prostate_output.log", "w")
# sys.stdout = logfile
# sys.stderr = logfile

# --- 2. Load and preprocess data ---
base_path = "../../ST_Datasets/Prostate"
ann_data = sc.read_visium(base_path)
print(ann_data)
data_module = make_spatial_data(ann_data)
data_module.prepare_data(rad_cutoff=300, n_top_genes=2000, min_cells=0)

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
model = STAMarker(5, "Prostate_STAMarker_output/", config)

# --- 4. Run the STAMarker pipeline ---
# Train autoencoders
model.train_auto_encoders(data_module)

# Perform clustering
n_clusters = 3
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
from sklearn.metrics import adjusted_rand_score
from scipy.spatial import cKDTree

# Load and prepare ground truth labels
all_true_labels_df = pd.read_csv("../../ST_Datasets/Prostate/prostate_true_labels.csv")
all_true_labels_df.rename(columns={'x': 'y', 'y': 'x'}, inplace=True) # Swap coordinates

# Align ground truth with experimental data using nearest neighbors
spatial_coords = pd.DataFrame(ann_data.obsm['spatial'], columns=['x', 'y'])
ari_true_labels_df = all_true_labels_df.dropna()
tree = cKDTree(spatial_coords[['x', 'y']])
distances, indices = tree.query(ari_true_labels_df[['x', 'y']], k=1)

# Calculate and print Adjusted Rand Index
consensus_aligned = ann_data.obs["Consensus clustering"].iloc[indices].to_numpy()
true_aligned = ari_true_labels_df['z'].to_numpy()
ari = adjusted_rand_score(true_aligned, consensus_aligned)
print(f"Adjusted Rand Index (ARI): {ari}")

# Generate and save comparison plot
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.embedding(ann_data, basis="spatial", s=25, color="Consensus clustering", ax=axes[0], show=False, title="Consensus Clustering")
temp_adata = ann_data[indices].copy()
temp_adata.obs['True Labels'] = ari_true_labels_df['z'].astype(str).values
sc.pl.embedding(temp_adata, basis="spatial", s=25, color="True Labels", ax=axes[1], show=False, title="True Labels")
plt.tight_layout()
plt.savefig("Prostate_STAMarker.png")

# --- 7. Generate final output CSV ---
# Create a dataframe with coordinates and consensus labels
consensus_df = pd.DataFrame(
    ann_data.obsm['spatial'],
    columns=['x', 'y'],
    index=ann_data.obs.index
)
consensus_df['barcode'] = consensus_df.index
consensus_df['STAMarker_label'] = ann_data.obs["Consensus clustering"].values

# Round coordinates for accurate merging
consensus_df['x'] = consensus_df['x'].round().astype(int)
consensus_df['y'] = consensus_df['y'].round().astype(int)
all_true_labels_df['x'] = all_true_labels_df['x'].round().astype(int)
all_true_labels_df['y'] = all_true_labels_df['y'].round().astype(int)

# Merge with ground truth labels
output_df = pd.merge(
    consensus_df,
    all_true_labels_df,
    on=['x', 'y'],
    how='left'
).rename(columns={'z': 'true_label'})

# Save to CSV
output_df.to_csv("STAMarker_Prostate_joint.csv", index=False)