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

csv_path = "/scratch/user/varogovchenko/BASTION_HPRC/STARmap/starmap_pcs_coords_labels.csv"
out_dir = "STARmap"
os.makedirs(out_dir, exist_ok=True)

logfile_path = os.path.join(out_dir, "SpaMask_STARmap_output.log")
print(f"Logging to {logfile_path}")

# Redirect stdout and stderr to the log file
original_stdout = sys.stdout
original_stderr = sys.stderr
logfile = open(logfile_path, "w")
sys.stdout = logfile
sys.stderr = logfile

print("=== SpaMask on STARmap (PC=10; k=7) ===")
print(f"Reading: {csv_path}")

# --------------------------- Load CSV ---------------------------------------
df = pd.read_csv(csv_path)

required_cols = ["ID", "x", "y"] + [f"PC{i}" for i in range(1, 11)]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in CSV: {missing}")

n = df.shape[0]
print(f"Rows (cells/spots): {n}")

# --------------------------- Build AnnData for SpaMask ------------------------
n_pcs = 10
pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
features = df[pc_cols].to_numpy(dtype=np.float32)

adata = ad.AnnData(np.zeros((n, 1), dtype=np.float32))
adata.obs_names = df["ID"].astype(str).values
adata.obsm["spatial"] = df[["x", "y"]].to_numpy(dtype=np.float32)
adata.obsm["feat"] = features

print(f"Feature matrix shape: {adata.obsm['feat'].shape}")

# --------------------------- Run SpaMask ------------------------------------
num_clusters = 7
k_cutoff = 12 

net = SPAMASK(adata,
             tissue_name='STARmap',
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

# If ground truth exists, compute ARI
if "ground_truth" in df.columns:
    gt = df["ground_truth"].astype('category').cat.codes
    pred = adata.obs[cluster_key].astype('category').cat.codes
    ari_val = metrics.adjusted_rand_score(gt, pred)
    print(f"ARI (n_pcs={n_pcs}, k={num_clusters}): {ari_val:.6f}")
    
    with open(os.path.join(out_dir, "starmap_spamask_ari.txt"), "w") as f:
        f.write(f"pca={n_pcs}, k={num_clusters} -> ARI={ari_val:.6f}\n")
    print("✓ Saved ARI summary to starmap_spamask_ari.txt")


csv_out = os.path.join(out_dir, f"SpaMask_STARmap_{n_pcs}PCs.csv")
out_df.to_csv(csv_out, index=False)
print(f"✓ Saved: {csv_out}")

# --------------------------- Plotting ---------------------------------------
if "ground_truth" in df.columns:
    plt.figure(figsize=(14, 6))

    # Plot predicted labels
    ax1 = plt.subplot(1, 2, 1)
    scatter1 = ax1.scatter(df["x"], df["y"], c=adata.obs[cluster_key].astype('category').cat.codes, cmap='viridis', s=5)
    ax1.set_title("SpaMask Predicted Labels")
    ax1.set_xlabel("X coordinate")
    ax1.set_ylabel("Y coordinate")
    ax1.set_aspect('equal', adjustable='box')
    ax1.legend(handles=scatter1.legend_elements()[0], labels=np.unique(adata.obs[cluster_key].astype('category')).tolist(), title="Clusters")

    # Plot ground truth labels
    ax2 = plt.subplot(1, 2, 2)
    scatter2 = ax2.scatter(df["x"], df["y"], c=df["ground_truth"].astype('category').cat.codes, cmap='viridis', s=5)
    ax2.set_title("Ground Truth Labels")
    ax2.set_xlabel("X coordinate")
    ax2.set_ylabel("Y coordinate")
    ax2.set_aspect('equal', adjustable='box')
    ax2.legend(handles=scatter2.legend_elements()[0], labels=np.unique(df["ground_truth"].astype('category')).tolist(), title="Clusters")


    plt.tight_layout()
    plot_path = os.path.join(out_dir, "spamask_starmap_clusters.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    # Use original stdout to print the plot path to the console
    print(f"✓ Saved plot to {plot_path}", file=original_stdout)


print("\n=== Done: SpaMask on STARmap for PC=10 with k=7 ===")

# Restore stdout and stderr
sys.stdout = original_stdout
sys.stderr = original_stderr
logfile.close()
