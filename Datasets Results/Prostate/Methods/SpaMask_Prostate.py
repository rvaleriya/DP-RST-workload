import os
import sys

# Add the SpaMask repo to the Python path to resolve the module not found error.
script_dir = os.path.dirname(os.path.realpath(__file__))
project_root = os.path.dirname(script_dir)
spamask_path = os.path.join(project_root, 'SpaMask')
sys.path.insert(0, spamask_path)

import warnings
import numpy as np
import pandas as pd
import torch
import anndata as ad
from SpaMask.spaMask import SPAMASK
from SpaMask.utils import fix_seed
from sklearn import metrics
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')

# --------------------------- Config -----------------------------------------
random_seed = 42
fix_seed(random_seed)
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Device: {device}")

out_dir = "Prostate"
os.makedirs(out_dir, exist_ok=True)

logfile_path = os.path.join(out_dir, "SpaMask_Prostate_output.log")
print(f"Logging to {logfile_path}")

# Redirect stdout and stderr to the log file
original_stdout = sys.stdout
original_stderr = sys.stderr
logfile = open(logfile_path, "w")
sys.stdout = logfile
sys.stderr = logfile

print("=== SpaMask on Prostate (PC=10; k=3) ===")
print(f"Reading data from CSV files in {out_dir}...")

# --------------------------- Load Converted CSV data ---------------------------------------
loc_df = pd.read_csv(os.path.join(out_dir, "prostate_loc.csv"))
pcs_df = pd.read_csv(os.path.join(out_dir, "prostate_pcs.csv"))
labels_df = pd.read_csv(os.path.join(out_dir, "prostate_labels.csv"))

n = loc_df.shape[0]
print(f"Rows (cells/spots): {n}")

# --------------------------- Build AnnData for SpaMask ------------------------
n_pcs = 10
features = pcs_df.iloc[:, :n_pcs].values

# Construct the master DataFrame
df = pd.DataFrame({
    'ID': [f'spot_{i}' for i in range(n)],
    'x': loc_df['x'].values,
    'y': loc_df['y'].values,
    'ground_truth': labels_df['ground_truth'].values
})

adata = ad.AnnData(np.zeros((n, 1), dtype=np.float32))
adata.obs_names = df["ID"].astype(str).values
adata.obsm["spatial"] = df[["x", "y"]].to_numpy(dtype=np.float32)
adata.obsm["feat"] = features

print(f"Feature matrix shape: {adata.obsm['feat'].shape}")

# --------------------------- Run SpaMask ------------------------------------
num_clusters = 3
k_cutoff = 12 

net = SPAMASK(adata,
             tissue_name='Prostate',
             num_clusters=num_clusters,
             device=device,
             max_epoch=1000,
             hidden_dim=512,
             latent_dim=256,
             lam=2.0,
             feat_mask_rate=0.5,
             edge_drop_rate=0.2,
             random_seed=random_seed,
             k_cutoff=k_cutoff,
             graph_model='KNN'
             )
net.train()

method = "kmeans"
net.process(method=method)
adata = net.get_adata()

cluster_key = 'kmeans'
print(f"✓ Clustering done -> obs['{cluster_key}']")

# --------------------------- Save Results -----------------------------------
out_df = pd.DataFrame({
    "ID": df["ID"].astype(str).values,
    "x": df["x"].values,
    "y": df["y"].values,
    cluster_key: adata.obs[cluster_key].astype(str).values
})

# Filter out NA values for ARI calculation
df_filtered = df.dropna(subset=['ground_truth'])
non_na_indices = df_filtered.index

gt = df_filtered["ground_truth"].astype('category').cat.codes
pred = adata.obs[cluster_key][non_na_indices].astype('category').cat.codes
ari_val = metrics.adjusted_rand_score(gt, pred)
print(f"ARI (n_pcs={n_pcs}, k={num_clusters}): {ari_val:.6f}")

with open(os.path.join(out_dir, "prostate_spamask_ari.txt"), "w") as f:
    f.write(f"pca={n_pcs}, k={num_clusters} -> ARI={ari_val:.6f}\n")
print("✓ Saved ARI summary to prostate_spamask_ari.txt")

csv_out = os.path.join(out_dir, f"SpaMask_Prostate_{n_pcs}PCs.csv")
out_df.to_csv(csv_out, index=False)
print(f"✓ Saved: {csv_out}")

# --------------------------- Plotting ---------------------------------------
plt.figure(figsize=(14, 6))

# Plot predicted labels (all points)
ax1 = plt.subplot(1, 2, 1)
scatter1 = ax1.scatter(df["x"], df["y"], c=adata.obs[cluster_key].astype('category').cat.codes, cmap='viridis', s=5)
ax1.set_title("SpaMask Predicted Labels")
ax1.set_xlabel("X coordinate")
ax1.set_ylabel("Y coordinate")
ax1.set_aspect('equal', adjustable='box')
ax1.legend(handles=scatter1.legend_elements()[0], labels=np.unique(adata.obs[cluster_key].astype('category')).tolist(), title="Clusters")

# Plot ground truth labels (only non-NA points)
ax2 = plt.subplot(1, 2, 2)
scatter2 = ax2.scatter(df_filtered["x"], df_filtered["y"], c=gt, cmap='viridis', s=5)
ax2.set_title("Ground Truth Labels (Non-NA)")
ax2.set_xlabel("X coordinate")
ax2.set_ylabel("Y coordinate")
ax2.set_aspect('equal', adjustable='box')
ax2.legend(handles=scatter2.legend_elements()[0], labels=np.unique(df_filtered["ground_truth"].astype('category')).tolist(), title="Clusters")

plt.tight_layout()
plot_path = os.path.join(out_dir, "spamask_prostate_clusters.png")
plt.savefig(plot_path, dpi=300)
plt.close()
# Use original stdout to print the plot path to the console
print(f"✓ Saved plot to {plot_path}", file=original_stdout)

print("\n=== Done: SpaMask on Prostate for PC=10 with k=3 ===")

# Restore stdout and stderr
sys.stdout = original_stdout
sys.stderr = original_stderr
logfile.close()
