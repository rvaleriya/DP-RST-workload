#!/usr/bin/env python3
import os
import sys

# Ensure R_LIBS includes user's Rlibs (modules handle R_HOME and PATH)
user_rlibs = "/scratch/user/varogovchenko/Rlibs"
if 'R_LIBS' in os.environ:
    if user_rlibs not in os.environ['R_LIBS']:
        os.environ['R_LIBS'] = f"{user_rlibs}:{os.environ['R_LIBS']}"
else:
    # If R_LIBS not set by modules, set it with user's Rlibs
    r_home = os.environ.get('R_HOME', '/sw/eb/sw/R/4.4.2-gfbf-2024a')
    os.environ['R_LIBS'] = f"{user_rlibs}:{r_home}/lib64/R/library"
import torch
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
import random

# Set random seeds for reproducibility
seed = 42
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
if torch.cuda.is_available():
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

from GraphST import GraphST
sys.path.append('../GraphST/GraphST')
from utils import clustering

# Configuration
dataset_name = "osmFISH"
num_clusters = 11
h5ad_path = "/scratch/user/varogovchenko/BASTION_HPRC/osmFISH/osmfish_20251122011024.h5ad"
output_csv = "GraphST_osmFISH.csv"
output_ari = "GraphST_osmFISH_ARI.txt"
output_plot = "GraphST_osmFISH_plot.png"

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# Load data
print(f"\nLoading {dataset_name} dataset...")
adata = sc.read_h5ad(h5ad_path)
adata.var_names_make_unique()
print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

# osmFISH data is already log-normalized, skip normalization
# Use all genes (no HVG filtering for osmFISH)

# Find ground truth column
ground_truth_col = None
for col in ['ground_truth', 'true_labels', 'label', 'cluster', 'cell_type']:
    if col in adata.obs.columns:
        ground_truth_col = col
        break

# Prepare results dataframe
results = pd.DataFrame({'barcode': adata.obs_names}, index=adata.obs_names)
if 'spatial' in adata.obsm:
    results[['x', 'y']] = adata.obsm['spatial']

# Run GraphST and cluster with 3PCs
print("\n=== Running GraphST with 3 PCs ===")
model_3pc = GraphST.GraphST(adata, device=device)
adata_3pc = model_3pc.train()
clustering(adata_3pc, n_clusters=num_clusters, radius=50, method='mclust', refinement=True, n_pcs=3)
results['mclust_3PCs'] = adata_3pc.obs['mclust'].astype(str).values

# Run GraphST and cluster with 10PCs
print("\n=== Running GraphST with 10 PCs ===")
model_10pc = GraphST.GraphST(adata, device=device)
adata_10pc = model_10pc.train()
clustering(adata_10pc, n_clusters=num_clusters, radius=50, method='mclust', refinement=True, n_pcs=10)
results['mclust_10PCs'] = adata_10pc.obs['mclust'].astype(str).values

# Save results
results.to_csv(output_csv, index=False)
print(f"\nSaved results to {output_csv}")

# Calculate ARI scores
ari_3pc = None
ari_10pc = None

if ground_truth_col:
    print(f"\nComputing ARI using '{ground_truth_col}' as ground truth...")
    valid_mask = ~adata.obs[ground_truth_col].isna()
    
    if valid_mask.sum() > 0:
        gt = pd.Categorical(adata.obs.loc[valid_mask, ground_truth_col].astype(str)).codes
        
        pr_3pc = pd.Categorical(results.loc[valid_mask, 'mclust_3PCs'].astype(str)).codes
        ari_3pc = metrics.adjusted_rand_score(gt, pr_3pc)
        
        pr_10pc = pd.Categorical(results.loc[valid_mask, 'mclust_10PCs'].astype(str)).codes
        ari_10pc = metrics.adjusted_rand_score(gt, pr_10pc)
        
        print(f"ARI (3PCs): {ari_3pc:.6f}")
        print(f"ARI (10PCs): {ari_10pc:.6f}")
        print(f"Calculated on {valid_mask.sum()}/{len(adata)} observations")
    else:
        print("No valid ground truth labels found")

# Save ARI results
with open(output_ari, 'w') as f:
    f.write(f"GraphST clustering results for {dataset_name}\n")
    f.write("="*50 + "\n")
    f.write(f"Number of clusters: {num_clusters}\n")
    if ground_truth_col:
        f.write(f"Ground truth column: {ground_truth_col}\n")
    f.write("="*50 + "\n")
    if ari_3pc is not None:
        f.write(f"ARI (3PCs): {ari_3pc:.6f}\n")
    if ari_10pc is not None:
        f.write(f"ARI (10PCs): {ari_10pc:.6f}\n")
    if ari_3pc is None and ari_10pc is None:
        f.write("ARI calculation skipped: No valid ground truth labels found.\n")

print(f"Saved ARI results to {output_ari}")

# Create plot
fig, axs = plt.subplots(1, 3, figsize=(18, 5))

if 'spatial' in adata.obsm:
    # Plot 3PCs results
    axs[0].scatter(adata.obsm['spatial'][:, 0], adata.obsm['spatial'][:, 1],
                   c=pd.Categorical(results['mclust_3PCs']).codes, cmap='tab10', s=1)
    axs[0].set_title(f'{dataset_name} - 3PCs')
    axs[0].set_xlabel('X')
    axs[0].set_ylabel('Y')
    axs[0].set_aspect('equal')
    
    # Plot 10PCs results
    axs[1].scatter(adata.obsm['spatial'][:, 0], adata.obsm['spatial'][:, 1],
                   c=pd.Categorical(results['mclust_10PCs']).codes, cmap='tab10', s=1)
    axs[1].set_title(f'{dataset_name} - 10PCs')
    axs[1].set_xlabel('X')
    axs[1].set_ylabel('Y')
    axs[1].set_aspect('equal')
    
    # Plot true labels
    if ground_truth_col:
        valid_mask = ~adata.obs[ground_truth_col].isna()
        if valid_mask.sum() > 0:
            axs[2].scatter(adata.obsm['spatial'][valid_mask, 0], 
                          adata.obsm['spatial'][valid_mask, 1],
                          c=pd.Categorical(adata.obs.loc[valid_mask, ground_truth_col].astype(str)).codes,
                          cmap='tab10', s=1)
            axs[2].set_title(f'{dataset_name} - True Labels')
        else:
            axs[2].text(0.5, 0.5, 'No valid\nground truth', ha='center', va='center', transform=axs[2].transAxes)
            axs[2].set_title(f'{dataset_name} - True Labels (N/A)')
    else:
        axs[2].text(0.5, 0.5, 'No ground truth\navailable', ha='center', va='center', transform=axs[2].transAxes)
        axs[2].set_title(f'{dataset_name} - True Labels (N/A)')
    axs[2].set_xlabel('X')
    axs[2].set_ylabel('Y')
    axs[2].set_aspect('equal')

plt.tight_layout()
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved plot to {output_plot}")

print(f"\nFinished processing {dataset_name}")

