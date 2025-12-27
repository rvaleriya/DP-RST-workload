print("--- Python script starting ---")

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
dataset_name = "Lung"
num_clusters = 6
n_pcs = 3
h5ad_path = "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/Data_VUILD96MF/Lung_xenium_processed.h5ad"
output_csv = "GraphST_Lung_3PCs.csv"
output_ari = "GraphST_Lung_3PCs_ARI.txt"
output_plot = "GraphST_Lung_3PCs_plot.png"

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# Load data
print(f"\nLoading {dataset_name} dataset...")
adata = sc.read_h5ad(h5ad_path)
adata.var_names_make_unique()
print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

# Handle preprocessed data
if 'logcounts' in adata.layers:
    adata.X = adata.layers['logcounts'].copy()
elif 'counts' in adata.layers:
    adata.X = adata.layers['counts'].copy()

# Set up spatial coordinates
if 'spatial' not in adata.obsm:
    if 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
        adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
        print("Using x_centroid and y_centroid for spatial coordinates")

# Preprocessing
if 'highly_variable' not in adata.var.columns:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()
print(f"After HVG filtering: {adata.n_vars} genes")

# Find ground truth column
ground_truth_col = None
for col in ['Annotation_Type', 'annotation', 'ground_truth', 'true_labels', 'label', 'cell_type']:
    if col in adata.obs.columns:
        ground_truth_col = col
        break

# Prepare results dataframe
results = pd.DataFrame({'barcode': adata.obs_names}, index=adata.obs_names)
if 'spatial' in adata.obsm:
    results[['x', 'y']] = adata.obsm['spatial']

# Run GraphST and cluster with 3PCs
print(f"\n=== Running GraphST with {n_pcs} PCs ===")
model = GraphST.GraphST(adata, device=device)
adata_result = model.train()
clustering(adata_result, n_clusters=num_clusters, radius=50, method='mclust', refinement=True, n_pcs=n_pcs)
results['mclust'] = adata_result.obs['mclust'].astype(str).values

# Save results
results.to_csv(output_csv, index=False)
print(f"\nSaved results to {output_csv}")

# Calculate ARI scores
ari_score = None

if ground_truth_col:
    print(f"\nComputing ARI using '{ground_truth_col}' as ground truth...")
    valid_mask = ~adata.obs[ground_truth_col].isna()
    
    if valid_mask.sum() > 0:
        gt = pd.Categorical(adata.obs.loc[valid_mask, ground_truth_col].astype(str)).codes
        
        pr = pd.Categorical(results.loc[valid_mask, 'mclust'].astype(str)).codes
        ari_score = metrics.adjusted_rand_score(gt, pr)
        
        print(f"ARI ({n_pcs}PCs): {ari_score:.6f}")
        print(f"Calculated on {valid_mask.sum()}/{len(adata)} observations")
    else:
        print("No valid ground truth labels found")

# Save ARI results
with open(output_ari, 'w') as f:
    f.write(f"GraphST clustering results for {dataset_name}\n")
    f.write("="*50 + "\n")
    f.write(f"Number of clusters: {num_clusters}\n")
    f.write(f"Number of PCs: {n_pcs}\n")
    if ground_truth_col:
        f.write(f"Ground truth column: {ground_truth_col}\n")
    f.write("="*50 + "\n")
    if ari_score is not None:
        f.write(f"ARI ({n_pcs}PCs): {ari_score:.6f}\n")
    else:
        f.write("ARI calculation skipped: No valid ground truth labels found.\n")

print(f"Saved ARI results to {output_ari}")

# Create plot
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

if 'spatial' in adata.obsm:
    # Plot clustering results
    axs[0].scatter(adata.obsm['spatial'][:, 0], adata.obsm['spatial'][:, 1],
                   c=pd.Categorical(results['mclust']).codes, cmap='tab10', s=1)
    axs[0].set_title(f'{dataset_name} - {n_pcs}PCs')
    axs[0].set_xlabel('X')
    axs[0].set_ylabel('Y')
    axs[0].set_aspect('equal')
    
    # Plot true labels
    if ground_truth_col:
        valid_mask = ~adata.obs[ground_truth_col].isna()
        if valid_mask.sum() > 0:
            axs[1].scatter(adata.obsm['spatial'][valid_mask, 0], 
                          adata.obsm['spatial'][valid_mask, 1],
                          c=pd.Categorical(adata.obs.loc[valid_mask, ground_truth_col].astype(str)).codes,
                          cmap='tab10', s=1)
            axs[1].set_title(f'{dataset_name} - True Labels')
        else:
            axs[1].text(0.5, 0.5, 'No valid\nground truth', ha='center', va='center', transform=axs[1].transAxes)
            axs[1].set_title(f'{dataset_name} - True Labels (N/A)')
    else:
        axs[1].text(0.5, 0.5, 'No ground truth\navailable', ha='center', va='center', transform=axs[1].transAxes)
        axs[1].set_title(f'{dataset_name} - True Labels (N/A)')
    axs[1].set_xlabel('X')
    axs[1].set_ylabel('Y')
    axs[1].set_aspect('equal')

plt.tight_layout()
plt.savefig(output_plot, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved plot to {output_plot}")

print(f"\nFinished processing {dataset_name} with {n_pcs} PCs")

