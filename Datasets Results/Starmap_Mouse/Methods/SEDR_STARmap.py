import os
os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R'
os.environ['R_LIBS_USER'] = '/scratch/user/varogovchenko/Rlibs'

import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import torch

from anndata import AnnData
from sklearn import metrics

import SEDR

# --------------------------- Config -----------------------------------------
random_seed = 42
SEDR.fix_seed(random_seed)
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

csv_path = "/scratch/user/varogovchenko/BASTION_HPRC/STARmap/starmap_pcs_coords_labels.csv"
out_dir  = "res_starmap_csv"
os.makedirs(out_dir, exist_ok=True)

logfile = open("SEDR_STARmap_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

print("=== SEDR on STARmap (PC=3 and PC=10; k=7) ===")
print(f"Device: {device}")
print(f"Reading: {csv_path}")

# --------------------------- Load CSV ---------------------------------------
df = pd.read_csv(csv_path)

required_cols = ["ID", "x", "y", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in CSV: {missing}")

n = df.shape[0]
print(f"Rows (cells/spots): {n}")

# --------------------------- Build AnnData for graph ------------------------
# SEDR's graph builder uses adata.obsm['spatial'] for coordinates
X_dummy = np.zeros((n, 1), dtype=np.float32)
adata = AnnData(X_dummy)
adata.obs_names = df["ID"].astype(str).values
adata.obsm["spatial"] = df[["x", "y"]].to_numpy(dtype=np.float32)

# Spatial KNN graph from coords
k_neighbors = 12
graph_dict = SEDR.graph_construction(adata, k_neighbors)
print("✓ Spatial graph constructed (k=12).")

# --------------------------- Helper: run one PC setup -----------------------
def run_sedr_for_pcs(n_pcs: int, k_clusters: int = 7):
    assert n_pcs in (3, 10), "This script is set up for 3 or 10 PCs."
    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    feat = df[pc_cols].to_numpy(dtype=np.float32)

    sedr_key = f"SEDR_{n_pcs}"
    cluster_key = f"mclust_pca{n_pcs}_{k_clusters}"

    print(f"\n--- Running SEDR for n_pcs={n_pcs}, k={k_clusters} ---")
    print(f"Feature matrix shape: {feat.shape}")

    sedr_net = SEDR.Sedr(feat, graph_dict, mode='clustering', device=device)
    sedr_net.train_with_dec(N=1)
    sedr_feat, _, _, _ = sedr_net.process()
    adata.obsm[sedr_key] = sedr_feat
    print(f"✓ Stored embedding in adata.obsm['{sedr_key}'] with shape {sedr_feat.shape}")

    SEDR.mclust_R(adata, k_clusters, use_rep=sedr_key, key_added=cluster_key)
    print(f"✓ Clustering done -> obs['{cluster_key}']")

    # Prepare tidy output
    out = pd.DataFrame({
        "ID": df["ID"].astype(str).values,
        "x": df["x"].values,
        "y": df["y"].values,
        cluster_key: adata.obs[cluster_key].astype(str).values
    })

    # Optional ARI if ground truth present
    ari_val = None
    if "ground_truth" in df.columns:
        gt = pd.Categorical(df["ground_truth"].astype(str).values).codes
        pred = pd.Categorical(adata.obs[cluster_key].astype(str).values).codes
        ari_val = metrics.adjusted_rand_score(gt, pred)
        print(f"ARI (n_pcs={n_pcs}, k={k_clusters}): {ari_val:.6f}")

    # Save per-PC CSV
    csv_out = os.path.join(out_dir, f"SEDR_STARmap_{n_pcs}PCs_k{k_clusters}.csv")
    out.to_csv(csv_out, index=False)
    print(f"✓ Saved: {csv_out}")

    return cluster_key, ari_val

# --------------------------- Run both PC settings ---------------------------
keys_and_ari = []
for pcs in (3, 10):
    key, ari = run_sedr_for_pcs(pcs, k_clusters=7)
    keys_and_ari.append((pcs, key, ari))

# --------------------------- Save combined labels --------------------------
combined = pd.DataFrame({
    "ID": df["ID"].astype(str).values,
    "x": df["x"].values,
    "y": df["y"].values,
})

for pcs, key, _ in keys_and_ari:
    combined[f"label_pca{pcs}_k7"] = adata.obs[key].astype(str).values

combined_out = os.path.join(out_dir, "SEDR_STARmap_3PCs_10PCs.csv")
combined.to_csv(combined_out, index=False)
print(f"\n✓ Saved combined labels: {combined_out}")

# Save ARIs (if computed)
ari_lines = []
for pcs, key, ari in keys_and_ari:
    if ari is not None:
        ari_lines.append(f"pca={pcs}, k=7 -> ARI={ari:.6f}")
if ari_lines:
    with open(os.path.join(out_dir, "starmap_sedr_ari.txt"), "w") as f:
        for line in ari_lines:
            f.write(line + "\n")
    print("✓ Saved ARI summary to starmap_sedr_ari.txt")

print("\n=== Done: SEDR on STARmap for PC=3 and PC=10 with k=7 ===")
logfile.close()