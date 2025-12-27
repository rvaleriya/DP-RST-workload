print("--- Python script starting ---")

import os
import sys
import torch
import scanpy as sc
import pandas as pd
from sklearn import metrics

from GraphST import GraphST
sys.path.append('../GraphST/GraphST')
from utils import clustering

# ----------------- Config -----------------
H5AD = "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/Data_VUILD96MF/Lung_xenium_processed.h5ad"
OUTPUT_CSV = "GraphST_Lung_Xenium_labels.csv"
OUTPUT_ARI_TXT = "graphst_lung_xenium_ARI.txt"
N_CLUSTERS = 6  # Fallback value, will be determined from data if possible
K = N_CLUSTERS  # for printing consistency

# Device
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Device: {device}")

# ----------------- Load Data -----------------
adata = sc.read_h5ad(H5AD)
adata.var_names_make_unique()
print(adata)

# Use pre-processed data from logcounts layer
if 'logcounts' in adata.layers:
    print("Using 'counts' layer as main data.")
    adata.X = adata.layers['counts'].copy()
    print(adata.X.head())
else:
    print("⚠️ 'logcounts' layer not found. Using adata.X.")

# Set up spatial coordinates
if 'x_centroid' in adata.obs.columns and 'y_centroid' in adata.obs.columns:
    print("Setting up spatial coordinates from 'x_centroid' and 'y_centroid'.")
    adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
else:
    print("⚠️ Spatial coordinates 'x_centroid' and 'y_centroid' not found in adata.obs.")

# ----------------- Preprocessing -----------------
# Data is already pre-processed (normalized and log1p)

# Use pre-computed highly variable genes
if 'highly_variable' in adata.var.columns:
    print("Using pre-computed highly variable genes.")
    adata = adata[:, adata.var.highly_variable].copy()
else:
    print("⚠️ 'highly_variable' not in adata.var.columns. Using top 2000 highly variable genes.")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()

# Determine number of clusters from ground truth
if "Annotation_Type" in adata.obs.columns:
    # Exclude NA values for clustering count
    n_clusters_auto = adata.obs["Annotation_Type"].dropna().nunique()
    if n_clusters_auto > 1:
        N_CLUSTERS = n_clusters_auto
        K = N_CLUSTERS
        print(f"Number of clusters (K) determined from 'Annotation_Type': {N_CLUSTERS}")
    else:
        print(f"⚠️ Only {n_clusters_auto} unique non-NA value(s) in 'Annotation_Type'. Using fallback K={K}")
else:
    print(f"⚠️ 'Annotation_Type' not found in adata.obs. Using fallback K={K}")


# ----------------- Prepare Results Table -----------------
results = pd.DataFrame({'barcode': adata.obs_names})
if "spatial" in adata.obsm:
    results[['x', 'y']] = adata.obsm['spatial']
else:
    print("⚠️  adata.obsm['spatial'] is missing; x/y will not be saved.")
    results['x'] = pd.NA
    results['y'] = pd.NA

# ----------------- Run GraphST and Clustering -----------------
print(f"\n=== Running GraphST ===")

# Train GraphST
model = GraphST.GraphST(adata, device=device)
adata_run = model.train()

# Clustering with mclust via rpy2
clustering(adata_run, n_clusters=N_CLUSTERS, radius=50, method='mclust', refinement=True)

# Store predicted labels in both 'results' and back into adata.obs
col_name_csv = f"mclust_{N_CLUSTERS}"
col_name_obs = "label"

labels = adata_run.obs["mclust"].astype(str).values
results[col_name_csv] = labels

# Ensure obs are aligned; GraphST preserves order, so direct assignment is fine
adata.obs[col_name_obs] = labels

# ----------------- Save Labels to CSV -----------------
results.to_csv(OUTPUT_CSV, index=False)
print(f"✓ Saved all cluster labels to {OUTPUT_CSV}")

# ----------------- Compute and Save ARI (if ground truth exists) -----------------
if "Annotation_Type" in adata.obs.columns:
    print("\nComputing ARI scores...")

    # Filter out NA values from ground truth for ARI calculation
    valid_cells = adata.obs['Annotation_Type'].notna()
    adata_filt = adata[valid_cells, :].copy()
    
    if adata_filt.n_obs == 0:
        print("⚠️ No valid ground truth labels found after filtering NAs. Skipping ARI.")
    else:
        gt = pd.Categorical(adata_filt.obs["Annotation_Type"].astype(str)).codes

        ari_lines = []
        print("\n" + "="*50)
        print("RESULTS (ARI on non-NA Annotation_Type):")
        print("="*50)

        if "label" in adata_filt.obs.columns:
            pr = pd.Categorical(adata_filt.obs["label"]).codes
            ari = metrics.adjusted_rand_score(gt, pr)
            print(f"  ARI (k={K}):  {ari:.4f}")
            ari_lines.append(f"ARI = {ari:.6f}")
        else:
            print("  ARI: missing label — skipped")

        print("="*50)

        # Save ARI summary
        with open(OUTPUT_ARI_TXT, "w") as f:
            f.write(f"GraphST + mclust clustering results (k={K}):\n")
            f.write("="*50 + "\n")
            if ari_lines:
                for line in ari_lines:
                    f.write(line + "\n")
            else:
                f.write("No ARI computed (missing labels or ground truth).\n")

        print(f"✅ Saved ARI summary: {OUTPUT_ARI_TXT}")
else:
    print("\n⚠️  No 'Annotation_Type' column found — ARI computation skipped.")

print("\nDone.")
