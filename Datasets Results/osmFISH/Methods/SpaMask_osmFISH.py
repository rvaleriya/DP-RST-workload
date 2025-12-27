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

csv_path = "/scratch/user/varogovchenko/BASTION_HPRC/osmFISH/osmfish_pcs_coords_labels.csv"
out_dir = "osmFISH"
os.makedirs(out_dir, exist_ok=True)

logfile_path = os.path.join(out_dir, "SpaMask_osmFISH_output.log")
print(f"Logging to {logfile_path}")

# Redirect stdout and stderr to the log file
original_stdout = sys.stdout
original_stderr = sys.stderr
logfile = open(logfile_path, "w")
sys.stdout = logfile
sys.stderr = logfile

print("=== SpaMask on osmFISH (running for 3PCs and 10PCs) ===")
print(f"Reading: {csv_path}")

# --------------------------- Load CSV ---------------------------------------
df = pd.read_csv(csv_path)

# Check for required columns (need up to PC10)
required_cols = ["ID", "x", "y"] + [f"PC{i}" for i in range(1, 11)]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in CSV: {missing}")

n = df.shape[0]
print(f"Rows (cells/spots): {n}")

# --------------------------- Run for both 3PCs and 10PCs ------------------------
num_clusters = 11
k_cutoff = 12

for n_pcs in [3, 10]:
    print(f"\n{'='*60}")
    print(f"=== Running SpaMask with {n_pcs} PCs ===")
    print(f"{'='*60}\n")
    
    # Build AnnData for SpaMask
    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    features = df[pc_cols].to_numpy(dtype=np.float32)
    
    adata = ad.AnnData(np.zeros((n, 1), dtype=np.float32))
    adata.obs_names = df["ID"].astype(str).values
    adata.obsm["spatial"] = df[["x", "y"]].to_numpy(dtype=np.float32)
    adata.obsm["feat"] = features
    
    print(f"Feature matrix shape: {adata.obsm['feat'].shape}")
    
    # Run SpaMask
    net = SPAMASK(adata,
                 tissue_name='osmFISH',
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
    
    # Save Results
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
        
        ari_file = os.path.join(out_dir, f"osmfish_spamask_{n_pcs}PCs_ari.txt")
        with open(ari_file, "w") as f:
            f.write(f"pca={n_pcs}, k={num_clusters} -> ARI={ari_val:.6f}\n")
        print(f"✓ Saved ARI summary to osmfish_spamask_{n_pcs}PCs_ari.txt")
    
    csv_out = os.path.join(out_dir, f"SpaMask_osmFISH_{n_pcs}PCs.csv")
    out_df.to_csv(csv_out, index=False)
    print(f"✓ Saved: {csv_out}")
    
    # Plotting
    if "ground_truth" in df.columns:
        plt.figure(figsize=(14, 6))
        
        # Plot predicted labels
        ax1 = plt.subplot(1, 2, 1)
        scatter1 = ax1.scatter(df["x"], df["y"], c=adata.obs[cluster_key].astype('category').cat.codes, cmap='viridis', s=5)
        ax1.set_title(f"SpaMask Predicted Labels ({n_pcs} PCs)")
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
        plot_path = os.path.join(out_dir, f"spamask_osmfish_{n_pcs}PCs_clusters.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()
        print(f"✓ Saved plot to {plot_path}")
    
    print(f"\n=== Done: SpaMask on osmFISH for PC={n_pcs} with k={num_clusters} ===\n")

print("\n" + "="*60)
print("=== All runs completed ===")
print("="*60)

# Restore stdout and stderr
sys.stdout = original_stdout
sys.stderr = original_stderr
logfile.close()

print(f"✓ All results saved to {out_dir}/", file=original_stdout)


