import os
import subprocess

# ============================
# R CONFIGURATION (must be set before any rpy2 imports)
# ============================
try:
    r_home_output = subprocess.check_output(['R', 'RHOME'], text=True, stderr=subprocess.DEVNULL).strip()
    if r_home_output and os.path.exists(r_home_output):
        os.environ['R_HOME'] = r_home_output
    else:
        os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R'
except (subprocess.CalledProcessError, FileNotFoundError):
    os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R'

os.environ['R_LIBS_USER'] = '/scratch/user/varogovchenko/Rlibs'

try:
    rscript_path = subprocess.check_output(['which', 'Rscript'], text=True).strip()
    if rscript_path:
        r_bin = os.path.dirname(rscript_path)
        if r_bin not in os.environ.get('PATH', ''):
            os.environ['PATH'] = r_bin + os.pathsep + os.environ.get('PATH', '')
except (subprocess.CalledProcessError, FileNotFoundError):
    r_bin = '/sw/eb/sw/R/4.4.2-gfbf-2024a/bin'
    if r_bin not in os.environ.get('PATH', ''):
        os.environ['PATH'] = r_bin + os.pathsep + os.environ.get('PATH', '')

import sys
from pathlib import Path
# Add parent directory to path to find SEDR module
sys.path.insert(0, str(Path(__file__).parent.parent))
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

dataset_name = "STARmap"
num_clusters = 7
csv_path = "/scratch/user/varogovchenko/BASTION_HPRC/STARmap/starmap_pcs_coords_labels.csv"

logfile = open("SEDR_STARmap_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

print(f"=== SEDR on {dataset_name} (PC=3 and PC=10; k={num_clusters}) ===")
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
X_dummy = np.zeros((n, 1), dtype=np.float32)
adata = AnnData(X_dummy)
adata.obs_names = df["ID"].astype(str).values
adata.obsm["spatial"] = df[["x", "y"]].to_numpy(dtype=np.float32)

# Spatial KNN graph from coords
k_neighbors = 12
graph_dict = SEDR.graph_construction(adata, k_neighbors)
print("✓ Spatial graph constructed (k=12).")

# --------------------------- Helper: run one PC setup -----------------------
def run_sedr_for_pcs(n_pcs: int, k_clusters: int = num_clusters):
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

    # Calculate ARI if ground truth is available
    ari_val = None
    if "ground_truth" in df.columns:
        try:
            gt = pd.Categorical(df["ground_truth"].astype(str).values).codes
            pred = pd.Categorical(adata.obs[cluster_key].astype(str).values).codes
            
            print(f"\nARI calculation (using 'ground_truth' as ground truth):")
            print(f"  Number of cells: {len(gt)}")
            print(f"  First 10 IDs: {list(df['ID'].astype(str)[:10])}")
            print(f"  First 10 ground truth labels: {list(df['ground_truth'].astype(str)[:10])}")
            print(f"  First 10 predicted clusters: {list(adata.obs[cluster_key].astype(str)[:10])}")
            print(f"  Unique ground truth values: {sorted(set(df['ground_truth'].astype(str)))}")
            print(f"  Unique predicted clusters: {sorted(set(adata.obs[cluster_key].astype(str)))}")
            
            # Create contingency table
            gt_labels = df["ground_truth"].astype(str).values
            pred_labels = adata.obs[cluster_key].astype(str).values
            contingency = pd.crosstab(
                pd.Series(gt_labels, name="Ground Truth"),
                pd.Series(pred_labels, name="Predicted")
            )
            print(f"\n  Contingency table (rows=GT, cols=Predicted):")
            print(contingency.to_string())
            
            ari_val = metrics.adjusted_rand_score(gt, pred)
            print(f"\nARI (n_pcs={n_pcs}, k={k_clusters}): {ari_val:.6f}")
        except Exception as e:
            print(f"Error calculating ARI: {e}")
            ari_val = None
    else:
        print("No 'ground_truth' column found in CSV. Skipping ARI calculation.")

    return cluster_key, ari_val

# --------------------------- Run both PC settings ---------------------------
cluster_keys_and_ari = []
for pcs in (3, 10):
    key, ari = run_sedr_for_pcs(pcs, k_clusters=num_clusters)
    cluster_keys_and_ari.append((pcs, key, ari))

# --------------------------- Save combined results --------------------------
print("\n--- Saving combined results ---")
combined = pd.DataFrame({
    "barcode": df["ID"].astype(str).values,
    "x": df["x"].values,
    "y": df["y"].values,
})

for pcs, key, _ in cluster_keys_and_ari:
    combined[f"label_pca{pcs}_k{num_clusters}"] = adata.obs[key].astype(str).values

# Add ground truth if available
if "ground_truth" in df.columns:
    combined['ground_truth'] = df["ground_truth"].astype(str).values

combined_out = "SEDR_STARmap.csv"
combined.to_csv(combined_out, index=False)
print(f"✓ Saved combined results: {combined_out}")

# Save ARI summary if computed
ari_lines = []
for pcs, key, ari in cluster_keys_and_ari:
    if ari is not None:
        ari_lines.append(f"pca={pcs}, k={num_clusters} -> ARI={ari:.6f}")
if ari_lines:
    print("\n--- ARI Summary ---")
    for line in ari_lines:
        print(line)
    with open("starmap_sedr_ari.txt", "w") as f:
        for line in ari_lines:
            f.write(line + "\n")
    print("✓ Saved ARI summary to starmap_sedr_ari.txt")

print(f"\n=== Done: SEDR on {dataset_name} for PC=3 and PC=10 with k={num_clusters} ===")
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__
logfile.close()

