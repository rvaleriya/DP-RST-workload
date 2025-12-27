import os
import subprocess

# ============================
# R CONFIGURATION (must be set before any rpy2 imports)
# ============================
# Dynamically detect R_HOME using R RHOME command (more reliable)
try:
    r_home_output = subprocess.check_output(['R', 'RHOME'], text=True, stderr=subprocess.DEVNULL).strip()
    if r_home_output and os.path.exists(r_home_output):
        os.environ['R_HOME'] = r_home_output
    else:
        # Fallback to known path
        os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R'
except (subprocess.CalledProcessError, FileNotFoundError):
    # Fallback to known path if R command not found
    os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R'

os.environ['R_LIBS_USER'] = '/scratch/user/varogovchenko/Rlibs'

# Ensure R bin directory is in PATH so rpy2 can find Rscript
try:
    rscript_path = subprocess.check_output(['which', 'Rscript'], text=True).strip()
    if rscript_path:
        r_bin = os.path.dirname(rscript_path)
        if r_bin not in os.environ.get('PATH', ''):
            os.environ['PATH'] = r_bin + os.pathsep + os.environ.get('PATH', '')
except (subprocess.CalledProcessError, FileNotFoundError):
    # Fallback to known path
    r_bin = '/sw/eb/sw/R/4.4.2-gfbf-2024a/bin'
    if r_bin not in os.environ.get('PATH', ''):
        os.environ['PATH'] = r_bin + os.pathsep + os.environ.get('PATH', '')

import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData
from sklearn import metrics

import SEDR

# --------------------------- Config -----------------------------------------
random_seed = 42
SEDR.fix_seed(random_seed)
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

csv_path = "/scratch/user/varogovchenko/BASTION_HPRC/osmFISH/osmfish_pcs_coords_labels.csv"
out_dir = "res_osmfish_csv"
os.makedirs(out_dir, exist_ok=True)

logfile = open("SEDR_osmFISH_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

print("=== SEDR on osmFISH (PC=3 and PC=10; k=11) ===")
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
def run_sedr_for_pcs(n_pcs: int, k_clusters: int = 11):
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
        
        # Debug: print first entries and contingency table
        print(f"\nDebug: ARI calculation check:")
        print(f"  Number of cells: {len(gt)}")
        print(f"  First 10 IDs: {list(df['ID'].astype(str)[:10])}")
        print(f"  First 10 ground truth labels: {list(df['ground_truth'].astype(str)[:10])}")
        print(f"  First 10 predicted clusters: {list(adata.obs[cluster_key].astype(str)[:10])}")
        print(f"  Unique ground truth values: {sorted(set(df['ground_truth'].astype(str)))}")
        print(f"  Unique predicted clusters: {sorted(set(adata.obs[cluster_key].astype(str)))}")
        
        # Create contingency table to see distribution
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

    # Save per-PC CSV
    csv_out = os.path.join(out_dir, f"SEDR_osmFISH_{n_pcs}PCs_k{k_clusters}.csv")
    out.to_csv(csv_out, index=False)
    print(f"✓ Saved: {csv_out}")

    # Create visualization if ground truth exists
    if "ground_truth" in df.columns:
        # Create AnnData for plotting (use dummy X)
        n_cells = len(df)
        plot_X = np.zeros((n_cells, 1), dtype=np.float32)
        plot_adata = AnnData(plot_X)
        plot_adata.obs_names = df["ID"].astype(str).values
        plot_adata.obsm["spatial"] = df[["x", "y"]].to_numpy(dtype=np.float32)
        plot_adata.obs[cluster_key] = adata.obs[cluster_key].astype(str).values
        plot_adata.obs["ground_truth"] = df["ground_truth"].astype(str).values
        
        # Create plot
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        sc.pl.embedding(plot_adata, basis="spatial", color=cluster_key, ax=axes[0], 
                       show=False, title=f"Predicted Clusters (PC={n_pcs}, k={k_clusters})", 
                       s=25, frameon=False)
        sc.pl.embedding(plot_adata, basis="spatial", color="ground_truth", ax=axes[1], 
                       show=False, title="Ground Truth Labels", s=25, frameon=False)
        plt.tight_layout()
        plot_out = os.path.join(out_dir, f"SEDR_osmFISH_{n_pcs}PCs_k{k_clusters}.png")
        plt.savefig(plot_out, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved plot: {plot_out}")

    return cluster_key, ari_val

# --------------------------- Run both PC settings ---------------------------
keys_and_ari = []
for pcs in (3, 10):
    key, ari = run_sedr_for_pcs(pcs, k_clusters=11)
    keys_and_ari.append((pcs, key, ari))

# --------------------------- Save combined labels --------------------------
combined = pd.DataFrame({
    "ID": df["ID"].astype(str).values,
    "x": df["x"].values,
    "y": df["y"].values,
})

for pcs, key, _ in keys_and_ari:
    combined[f"label_pca{pcs}_k11"] = adata.obs[key].astype(str).values

if "ground_truth" in df.columns:
    combined['ground_truth'] = df["ground_truth"].astype(str).values

combined_out = os.path.join(out_dir, "SEDR_osmFISH_3PCs_10PCs.csv")
combined.to_csv(combined_out, index=False)
print(f"\n✓ Saved combined labels: {combined_out}")

# Save ARIs (if computed)
ari_lines = []
for pcs, key, ari in keys_and_ari:
    if ari is not None:
        ari_lines.append(f"pca={pcs}, k=11 -> ARI={ari:.6f}")
if ari_lines:
    with open(os.path.join(out_dir, "osmfish_sedr_ari.txt"), "w") as f:
        for line in ari_lines:
            f.write(line + "\n")
    print("✓ Saved ARI summary to osmfish_sedr_ari.txt")

print("\n=== Done: SEDR on osmFISH for PC=3 and PC=10 with k=11 ===")
logfile.close()

