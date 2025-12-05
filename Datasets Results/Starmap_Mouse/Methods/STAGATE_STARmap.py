#!/usr/bin/env python3
import sys
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'  # Disable GPU
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'  # Reduce TensorFlow logging

# Add STAGATE_local_code to Python path
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
stagate_path = os.path.join(project_root, 'STAGATE_local_code')
sys.path.insert(0, stagate_path)

# Set CUDA paths
os.environ['CUDA_HOME'] = '/sw/eb/sw/CUDA/11.7.0'
os.environ['PATH'] = f"{os.environ['CUDA_HOME']}/bin:{os.environ['PATH']}"
os.environ['LD_LIBRARY_PATH'] = f"{os.environ['CUDA_HOME']}/lib64:{os.environ.get('LD_LIBRARY_PATH', '')}"

# Set R paths
os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a'
os.environ['R_LIBS'] = f"/scratch/user/varogovchenko/Rlibs:/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R/library"
os.environ['PATH'] = f"{os.environ['R_HOME']}/bin:{os.environ['PATH']}"

# Open log file in the same directory as the script
logfile = open("STAGATE_STARmap_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import torch
import tensorflow as tf
from sklearn import metrics
import random
import scipy.sparse as sp

# Set seeds for reproducibility
seed = 42
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
tf.random.set_seed(seed)
if torch.cuda.is_available():
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

tf.compat.v1.disable_eager_execution()
physical_devices = tf.config.list_physical_devices('GPU')
if physical_devices:
    try:
        for device in physical_devices:
            tf.config.experimental.set_memory_growth(device, True)
        print(f"Found {len(physical_devices)} GPU(s)")
    except RuntimeError as e:
        print(f"Error configuring GPU: {e}")
else:
    print("No GPU devices found")

from STAGATE.utils import Cal_Spatial_Net, Stats_Spatial_Net, mclust_R
from STAGATE.Train_STAGATE import train_STAGATE

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Device: {device}")

# Dataset configuration
dataset_name = "STARmap"
num_clusters = 7  # True number of clusters
h5ad_path = "/scratch/user/varogovchenko/BASTION_HPRC/STARmap/STARmap_20180505_BY3_1k_20251008011714.h5ad"

print(f"\n=== Processing {dataset_name} ===")
print(f"H5AD path: {h5ad_path}")

# Load data
adata = sc.read_h5ad(h5ad_path)
adata.var_names_make_unique()
print(adata)

# Preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()
print(adata)

# Convert to sparse matrix if needed (STAGATE requires sparse format)
if not sp.issparse(adata.X):
    print("Converting adata.X to sparse CSR format for STAGATE...")
    adata.X = sp.csr_matrix(adata.X)

print("Running STAGATE...")
Cal_Spatial_Net(adata, rad_cutoff=150)
Stats_Spatial_Net(adata)
adata = train_STAGATE(adata, alpha=0, hidden_dims=[512, 30])
sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)
adata = mclust_R(adata, used_obsm='STAGATE', num_cluster=num_clusters)
adata.obs['mclust'] = adata.obs['mclust'].astype(str)

# Save clustering results
df = pd.DataFrame({
    'barcode': adata.obs_names,
    'mclust': adata.obs['mclust'].values
})
if 'spatial' in adata.obsm:
    df[['x', 'y']] = adata.obsm['spatial']
else:
    print("Warning: No spatial coordinates found in obsm['spatial']")
    df['x'] = pd.NA
    df['y'] = pd.NA

# Save results
base_dir = os.path.dirname(os.path.abspath(__file__))
results_file = os.path.join(base_dir, "STAGATE_STARmap.csv")
df.to_csv(results_file, index=False)
print(f"\n✓ Saved clustering results to: {results_file}")

# Compute ARI scores
ari = None
ground_truth_col = None

# Check for ground truth column (common names)
for gt_col in ['ground_truth', 'true_labels', 'label', 'cluster', 'cell_type']:
    if gt_col in adata.obs.columns:
        ground_truth_col = gt_col
        break

if ground_truth_col:
    print(f"\n=== Computing ARI scores (using '{ground_truth_col}' as ground truth) ===")
    # Filter out NA values before calculating ARI
    valid_mask = ~adata.obs[ground_truth_col].isna()
    if valid_mask.sum() > 0:
        gt = pd.Categorical(adata.obs.loc[valid_mask, ground_truth_col].astype(str)).codes
        pr = pd.Categorical(adata.obs.loc[valid_mask, 'mclust']).codes
        
        ari = metrics.adjusted_rand_score(gt, pr)
        print(f"ARI: {ari:.6f}")
        if valid_mask.sum() < len(adata):
            print(f"Note: Calculated ARI on {valid_mask.sum()} out of {len(adata)} observations (excluded {len(adata) - valid_mask.sum()} NA values)")
    else:
        print("⚠️  No valid ground truth labels found (all are NA) — ARI calculation skipped")
    
    # Save ARI to text file
    ari_file = os.path.join(base_dir, "STAGATE_STARmap_ARI.txt")
    with open(ari_file, "w") as f:
        f.write(f"STAGATE clustering results for {dataset_name}\n")
        f.write("="*50 + "\n")
        f.write(f"Ground truth column: {ground_truth_col}\n")
        f.write(f"Number of clusters: {num_clusters}\n")
        f.write("="*50 + "\n")
        if ari is not None:
            f.write(f"ARI: {ari:.6f}\n")
        else:
            f.write("ARI calculation skipped: No valid ground truth labels found.\n")
    print(f"✓ Saved ARI results to: {ari_file}")
else:
    print("\n⚠️  No ground truth column found — ARI computation skipped.")
    # Create empty ARI file
    ari_file = os.path.join(base_dir, "STAGATE_STARmap_ARI.txt")
    with open(ari_file, "w") as f:
        f.write(f"STAGATE clustering results for {dataset_name}\n")
        f.write("="*50 + "\n")
        f.write("No ground truth column found — ARI computation skipped.\n")
    print(f"✓ Created ARI file (no ground truth): {ari_file}")

if 'spatial' in adata.obsm:
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    
    spatial_coords = adata.obsm['spatial']
    has_visium_spatial = 'spatial' in adata.uns and isinstance(adata.uns['spatial'], dict)
    
    if has_visium_spatial:
        sc.pl.spatial(adata, color='mclust', ax=axs[0], show=False,
                     title=f'{dataset_name} - STAGATE', spot_size=1.5)
        
        if ground_truth_col:
            sc.pl.spatial(adata, color=ground_truth_col, ax=axs[1], show=False,
                         title=f'{dataset_name} - True Labels', spot_size=1.5)
        else:
            axs[1].text(0.5, 0.5, 'No ground truth\navailable', 
                       ha='center', va='center', transform=axs[1].transAxes)
            axs[1].set_title(f'{dataset_name} - True Labels (N/A)')
    else:
        scatter1 = axs[0].scatter(spatial_coords[:, 0], spatial_coords[:, 1], 
                                  c=pd.Categorical(adata.obs['mclust']).codes, 
                                  cmap='viridis', s=5)
        axs[0].set_title(f'{dataset_name} - STAGATE')
        axs[0].set_xlabel('X coordinate')
        axs[0].set_ylabel('Y coordinate')
        axs[0].set_aspect('equal', adjustable='box')
        
        if ground_truth_col:
            scatter2 = axs[1].scatter(spatial_coords[:, 0], spatial_coords[:, 1], 
                                      c=pd.Categorical(adata.obs[ground_truth_col]).codes, 
                                      cmap='viridis', s=5)
            axs[1].set_title(f'{dataset_name} - True Labels')
            axs[1].set_xlabel('X coordinate')
            axs[1].set_ylabel('Y coordinate')
            axs[1].set_aspect('equal', adjustable='box')
        else:
            axs[1].text(0.5, 0.5, 'No ground truth\navailable', 
                       ha='center', va='center', transform=axs[1].transAxes)
            axs[1].set_title(f'{dataset_name} - True Labels (N/A)')
    
    plt.tight_layout()
    plot_file = os.path.join(base_dir, "STAGATE_STARmap_plot.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to: {plot_file}")
else:
    print("No spatial coordinates found — skipping spatial plot.")

print(f"\nFinished {dataset_name}")

# Close the log file at the end
logfile.close()
