#!/usr/bin/env python3
import sys
import os
import random
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
import SpaGCN as spg
from scipy.sparse import issparse
import torch
from pathlib import Path
from sklearn.metrics import adjusted_rand_score as ari
import matplotlib.pyplot as plt

# =========================
# Logging
# =========================
LOGFILE = "SpaGCN_osmFISH_output.log"
logfile = open(LOGFILE, "w")
sys.stdout = logfile
sys.stderr = logfile
print(f"--- Script Start: Logging to {LOGFILE} ---")

# =========================
# Config
# =========================
H5AD = "/scratch/user/varogovchenko/BASTION_HPRC/osmFISH/osmfish_20251122011024.h5ad"
OUTPUT_CSV = "SpaGCN_osmFISH.csv"
OUTPUT_ARI_TXT = "SpaGCN_osmFISH_ARI.txt"
OUTPUT_PLOT = "SpaGCN_osmFISH_plot.png"

N_CLUSTERS = 11
NUM_PCS_LIST = [3, 10]
P_PARAM = 0.5
SEED = 42

# =========================
# Environment / seeds
# =========================
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print("SpaGCN version:", getattr(spg, "__version__", "unknown"))
print("Device:", device)

# =========================
# Load osmFISH data
# =========================
adata = sc.read_h5ad(H5AD)
adata.var_names_make_unique()
print(adata)

# Coordinates
if "spatial" not in adata.obsm:
    raise RuntimeError("osmFISH AnnData lacks obsm['spatial'] coordinates.")
coords = adata.obsm["spatial"]
if coords.shape[1] < 2:
    raise RuntimeError("obsm['spatial'] must have at least 2 columns (x, y).")

adata.obs["x_array"] = coords[:, 0]
adata.obs["y_array"] = coords[:, 1]
adata.obs["x_pixel"] = coords[:, 0]
adata.obs["y_pixel"] = coords[:, 1]

# =========================
# Preprocessing (osmFISH data is already log-normalized)
# =========================
# Skip normalization and log transformation, keep all genes
# But still do PCA
max_pcs = max(NUM_PCS_LIST)
sc.tl.pca(adata, n_comps=max_pcs)
print(f"✓ Computed {max_pcs} PCs")

# =========================
# Adjacency (no histology)
# =========================
print("Calculating adjacency (histology=False)...")
x_array = adata.obs["x_array"].tolist()
y_array = adata.obs["y_array"].tolist()
adj = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)

# =========================
# Search l for p
# =========================
print(f"Searching l for p={P_PARAM} ...")
l = spg.search_l(P_PARAM, adj)
print("Recommended l =", l)

# =========================
# Check for ground truth
# =========================
ground_truth_col = None
for gt_col in ['ground_truth', 'true_labels', 'label', 'cluster', 'cell_type']:
    if gt_col in adata.obs.columns:
        ground_truth_col = gt_col
        break

# =========================
# Results container
# =========================
results = pd.DataFrame({
    "barcode": adata.obs_names,
    "x": coords[:, 0],
    "y": coords[:, 1],
})
ari_rows = []

# =========================
# Main loop over PC configs
# =========================
for num_pcs in NUM_PCS_LIST:
    print("\n" + "="*60)
    print(f"Running SpaGCN with {num_pcs} PCs (k={N_CLUSTERS})")
    print("="*60)

    adata_run = adata.copy()
    adata_run.obsm["X_pca"] = adata_run.obsm["X_pca"][:, :num_pcs]

    r_seed = t_seed = n_seed = SEED
    print(f"Using kmeans initialization (n_clusters={N_CLUSTERS}) to avoid neighbor computation issues...")
    
    # For osmFISH, use kmeans initialization to avoid neighbor computation issues
    # Kmeans doesn't require neighbor computation, so it avoids the broadcast error
    clf = spg.SpaGCN()
    clf.set_l(l)
    random.seed(r_seed); torch.manual_seed(t_seed); np.random.seed(n_seed)
    print("Training SpaGCN...")
    clf.train(
        adata_run, adj,
        init_spa=True, init="kmeans", n_clusters=N_CLUSTERS,
        tol=5e-3, lr=0.05, max_epochs=200, num_pcs=num_pcs
    )
    y_pred, prob = clf.predict()
    adata_run.obs["pred"] = pd.Categorical(y_pred)
    print("Initial prediction stored in adata_run.obs['pred']")

    print("Refining clusters (shape='square')...")
    adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    refined = spg.refine(
        sample_id=adata_run.obs.index.tolist(),
        pred=adata_run.obs["pred"].tolist(),
        dis=adj_2d,
        shape="square"
    )
    adata_run.obs["refined_pred"] = pd.Categorical(refined)
    print("Refined prediction stored in adata_run.obs['refined_pred']")

    col_init = f"pred_{num_pcs}PCs_k{N_CLUSTERS}"
    col_ref  = f"refined_pred_{num_pcs}PCs_k{N_CLUSTERS}"
    results[col_init] = adata_run.obs["pred"].astype(str).values
    results[col_ref]  = adata_run.obs["refined_pred"].astype(str).values

    # ARI calculation
    if ground_truth_col:
        valid_mask = ~adata_run.obs[ground_truth_col].isna()
        if valid_mask.sum() > 0:
            gt = adata_run.obs.loc[valid_mask, ground_truth_col].astype(str).tolist()
            pr_i = adata_run.obs.loc[valid_mask, "pred"].astype(str).tolist()
            pr_r = adata_run.obs.loc[valid_mask, "refined_pred"].astype(str).tolist()
            try:
                ari_i = ari(gt, pr_i)
            except Exception:
                ari_i = float("nan")
            try:
                ari_r = ari(gt, pr_r)
            except Exception:
                ari_r = float("nan")
            ari_rows.append({"PCs": num_pcs, "k": N_CLUSTERS, "ARI_initial": ari_i, "ARI_refined": ari_r})
            print(f"✓ ARI initial ({num_pcs} PCs, k={N_CLUSTERS}): {ari_i:.4f}")
            print(f"✓ ARI refined ({num_pcs} PCs, k={N_CLUSTERS}): {ari_r:.4f}")

# =========================
# Save results
# =========================
results.to_csv(OUTPUT_CSV, index=False)
print(f"\n✓ Saved clustering labels to {OUTPUT_CSV}")
print("Columns:", ", ".join(results.columns))

# Save ARI summary
if ari_rows:
    ari_df = pd.DataFrame(ari_rows)
    with open(OUTPUT_ARI_TXT, "w") as f:
        f.write("SpaGCN on osmFISH (k=11) — ARI summary\n")
        f.write("="*60 + "\n")
        for _, r in ari_df.iterrows():
            f.write(f"PCs: {int(r['PCs'])}  |  k: {int(r['k'])}  |  ARI initial: {r['ARI_initial']:.6f}  |  ARI refined: {r['ARI_refined']:.6f}\n")
        best_row = ari_df.loc[ari_df["ARI_refined"].idxmax()]
        f.write("\n" + "-"*60 + "\n")
        f.write(f"Best refined ARI: {best_row['ARI_refined']:.6f}  (PCs={int(best_row['PCs'])})\n")
    print(f"✓ Saved ARI summary to {OUTPUT_ARI_TXT}")
else:
    print("⚠️  No ground truth labels found — ARI summary not written.")

# =========================
# Create plot
# =========================
print("\nCreating visualization plot...")
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

spatial_coords = adata.obsm['spatial']

# 3PCs results
scatter1 = axes[0].scatter(spatial_coords[:, 0], spatial_coords[:, 1], 
                          c=pd.Categorical(results["refined_pred_3PCs_k11"]).codes, 
                          cmap='viridis', s=5)
axes[0].set_title(f'osmFISH - 3PCs (k={N_CLUSTERS})')
axes[0].set_xlabel('X coordinate')
axes[0].set_ylabel('Y coordinate')
axes[0].set_aspect('equal', adjustable='box')

# 10PCs results
scatter2 = axes[1].scatter(spatial_coords[:, 0], spatial_coords[:, 1], 
                          c=pd.Categorical(results["refined_pred_10PCs_k11"]).codes, 
                          cmap='viridis', s=5)
axes[1].set_title(f'osmFISH - 10PCs (k={N_CLUSTERS})')
axes[1].set_xlabel('X coordinate')
axes[1].set_ylabel('Y coordinate')
axes[1].set_aspect('equal', adjustable='box')

# True labels
if ground_truth_col and ground_truth_col in adata.obs.columns:
    valid_mask = ~adata.obs[ground_truth_col].isna()
    if valid_mask.sum() > 0:
        temp_coords = spatial_coords[valid_mask]
        temp_labels = adata.obs.loc[valid_mask, ground_truth_col]
        scatter3 = axes[2].scatter(temp_coords[:, 0], temp_coords[:, 1], 
                                  c=pd.Categorical(temp_labels).codes, 
                                  cmap='viridis', s=5)
        axes[2].set_title('osmFISH - True Labels')
    else:
        axes[2].text(0.5, 0.5, 'No valid ground truth\nlabels available', 
                   ha='center', va='center', transform=axes[2].transAxes)
        axes[2].set_title('osmFISH - True Labels (N/A)')
else:
    axes[2].text(0.5, 0.5, 'No ground truth\navailable', 
               ha='center', va='center', transform=axes[2].transAxes)
    axes[2].set_title('osmFISH - True Labels (N/A)')
axes[2].set_xlabel('X coordinate')
axes[2].set_ylabel('Y coordinate')
axes[2].set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved plot to {OUTPUT_PLOT}")

print("\nDone.")
logfile.close()

