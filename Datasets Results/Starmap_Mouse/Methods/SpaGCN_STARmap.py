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

# =========================
# Logging
# =========================
LOGFILE = "SpaGCN_STARmap_output.log"
logfile = open(LOGFILE, "w")
sys.stdout = logfile
sys.stderr = logfile
print(f"--- Script Start: Logging to {LOGFILE} ---")

# =========================
# Config (STARmap)
# =========================
H5AD = "/scratch/user/varogovchenko/BASTION_HPRC/STARmap/STARmap_20180505_BY3_1k_20251008011714.h5ad"
OUTPUT_CSV = "SpaGCN_STARmap_3PCs_10PCs.csv"
OUTPUT_ARI_TXT = "SpaGCN_STARmap_ARI.txt"

N_CLUSTERS = 7                    # STARmap: 7 clusters
NUM_PCS_LIST = [3, 10]            # Evaluate both PC settings
P_PARAM = 0.5                     # p for search_l
SEED = 42                         # reproducibility seed

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
# Load STARmap data
# =========================
adata = sc.read_h5ad(H5AD)
adata.var_names_make_unique()
print(adata)

# Coordinates (STARmap): expect obsm['spatial'] with (x, y)
if "spatial" not in adata.obsm:
    raise RuntimeError("STARmap AnnData lacks obsm['spatial'] coordinates.")
coords = adata.obsm["spatial"]
if coords.shape[1] < 2:
    raise RuntimeError("obsm['spatial'] must have at least 2 columns (x, y).")

# Map to SpaGCN-expected fields in .obs
adata.obs["x_array"] = coords[:, 0]
adata.obs["y_array"] = coords[:, 1]
adata.obs["x_pixel"] = coords[:, 0]   # No histology; duplicate for compatibility
adata.obs["y_pixel"] = coords[:, 1]

# =========================
# Preprocessing (SpaGCN style)
# =========================
# spg.prefilter_genes(adata, min_cells=3)   # if available in your SpaGCN
# spg.prefilter_specialgenes(adata)         # if available
# Normalize & log
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

# HVGs + PCA once
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
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

    # Prepare PCs for SpaGCN (SpaGCN reads adata.obsm['X_pca'])
    adata_run = adata.copy()
    adata_run.obsm["X_pca"] = adata_run.obsm["X_pca"][:, :num_pcs]

    # Search resolution for desired clusters
    r_seed = t_seed = n_seed = SEED
    print(f"Searching resolution for n_clusters={N_CLUSTERS} ...")
    # search_res will use adata.obsm['X_pca'] if present
    res = spg.search_res(
        adata_run, adj, l, N_CLUSTERS,
        start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20,
        r_seed=r_seed, t_seed=t_seed, n_seed=n_seed, num_pcs=num_pcs
    )
    print("Recommended res =", res)

    # Train SpaGCN
    clf = spg.SpaGCN()
    clf.set_l(l)
    random.seed(r_seed); torch.manual_seed(t_seed); np.random.seed(n_seed)
    print("Training SpaGCN...")
    clf.train(
        adata_run, adj,
        init_spa=True, init="louvain", res=res,
        tol=5e-3, lr=0.05, max_epochs=200, num_pcs=num_pcs
    )
    y_pred, prob = clf.predict()
    adata_run.obs["pred"] = pd.Categorical(y_pred)
    print("Initial prediction stored in adata_run.obs['pred']")

    # Refine (STARmap grid → 'square')
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

    # Save columns
    col_init = f"pred_{num_pcs}PCs_k{N_CLUSTERS}"
    col_ref  = f"refined_pred_{num_pcs}PCs_k{N_CLUSTERS}"
    results[col_init] = adata_run.obs["pred"].astype(str).values
    results[col_ref]  = adata_run.obs["refined_pred"].astype(str).values

    # ARI if ground_truth exists
    if "ground_truth" in adata_run.obs.columns:
        gt = adata_run.obs["ground_truth"].astype(str).tolist()
        pr_i = adata_run.obs["pred"].astype(str).tolist()
        pr_r = adata_run.obs["refined_pred"].astype(str).tolist()
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

# Save ARI summary (if any)
if ari_rows:
    ari_df = pd.DataFrame(ari_rows)
    with open(OUTPUT_ARI_TXT, "w") as f:
        f.write("SpaGCN on STARmap (k=7) — ARI summary\n")
        f.write("="*60 + "\n")
        for _, r in ari_df.iterrows():
            f.write(f"PCs: {int(r['PCs'])}  |  k: {int(r['k'])}  |  ARI initial: {r['ARI_initial']:.6f}  |  ARI refined: {r['ARI_refined']:.6f}\n")
        best_row = ari_df.loc[ari_df["ARI_refined"].idxmax()]
        f.write("\n" + "-"*60 + "\n")
        f.write(f"Best refined ARI: {best_row['ARI_refined']:.6f}  (PCs={int(best_row['PCs'])})\n")
    print(f"✓ Saved ARI summary to {OUTPUT_ARI_TXT}")
else:
    print("⚠️  No 'ground_truth' column found — ARI summary not written.")

print("\nDone.")
logfile.close()