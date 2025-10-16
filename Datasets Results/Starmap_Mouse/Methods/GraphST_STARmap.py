#!/usr/bin/env python3
import os
import sys
import torch
import scanpy as sc
import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt  # not used for plots now, but harmless if imported

from GraphST import GraphST
sys.path.append('GraphST/GraphST')
from utils import clustering

# ----------------- Config -----------------
H5AD = "/scratch/user/varogovchenko/BASTION_HPRC/STARmap/STARmap_20180505_BY3_1k_20251008011714.h5ad"
OUTPUT_CSV = "GraphST_STARmap_3PCs_10PCs.csv"
OUTPUT_ARI_TXT = "graphst_starmap_ARI.txt"
PCA_SETTINGS = [3, 10]
N_CLUSTERS = 7
K = N_CLUSTERS  # for printing consistency

# Redirect stdout/stderr to log file
logfile = open("GraphST_STARmap_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

# Device
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Device: {device}")

# ----------------- Load Data -----------------
adata = sc.read_h5ad(H5AD)
adata.var_names_make_unique()
print(adata)

# ----------------- Preprocessing -----------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()

# ----------------- Prepare Results Table -----------------
results = pd.DataFrame({'barcode': adata.obs_names})
if "spatial" in adata.obsm:
    results[['x', 'y']] = adata.obsm['spatial']
else:
    print("⚠️  adata.obsm['spatial'] is missing; x/y will not be saved.")
    results['x'] = pd.NA
    results['y'] = pd.NA

# ----------------- Run GraphST and Clustering for Each PC Setting -----------------
for n_pc in PCA_SETTINGS:
    print(f"\n=== Running GraphST with {n_pc} PCs ===")
    sc.pp.pca(adata, n_comps=n_pc)

    # Train GraphST
    model = GraphST.GraphST(adata, device=device)
    adata_run = model.train()

    # Clustering with mclust via rpy2
    clustering(adata_run, n_clusters=N_CLUSTERS, radius=50, method='mclust', refinement=True)

    # Store predicted labels in both 'results' and back into adata.obs
    col_name_csv = f"mclust_{N_CLUSTERS}_{n_pc}PCs"
    col_name_obs = f"label_{n_pc}PCs"

    labels = adata_run.obs["mclust"].astype(str).values
    results[col_name_csv] = labels

    # Ensure obs are aligned; GraphST preserves order, so direct assignment is fine
    adata.obs[col_name_obs] = labels

# ----------------- Save Labels to CSV -----------------
results.to_csv(OUTPUT_CSV, index=False)
print(f"✓ Saved all cluster labels to {OUTPUT_CSV}")

# ----------------- Compute and Save ARI (if ground truth exists) -----------------
if "ground_truth" in adata.obs.columns:
    print("\nComputing ARI scores...")

    gt = pd.Categorical(adata.obs["ground_truth"].astype(str)).codes

    ari_lines = []
    print("\n" + "="*50)
    print("RESULTS:")
    print("="*50)

    if "label_3PCs" in adata.obs.columns:
        pr3 = pd.Categorical(adata.obs["label_3PCs"]).codes
        ari3 = metrics.adjusted_rand_score(gt, pr3)
        print(f"  ARI (3 PCs, k={K}):  {ari3:.4f}")
        ari_lines.append(f"3PCs:  ARI = {ari3:.6f}")
    else:
        print("  ARI (3 PCs): missing label_3PCs — skipped")

    if "label_10PCs" in adata.obs.columns:
        pr10 = pd.Categorical(adata.obs["label_10PCs"]).codes
        ari10 = metrics.adjusted_rand_score(gt, pr10)
        print(f"  ARI (10 PCs, k={K}): {ari10:.4f}")
        ari_lines.append(f"10PCs: ARI = {ari10:.6f}")
    else:
        print("  ARI (10 PCs): missing label_10PCs — skipped")

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
    print("\n⚠️  No ground_truth column found — ARI computation skipped.")

print("\nDone.")
logfile.close()