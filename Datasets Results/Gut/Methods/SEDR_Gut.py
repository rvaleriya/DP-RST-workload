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
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn import metrics
from scipy.spatial import cKDTree

import SEDR

# --------------------------- Config -----------------------------------------
random_seed = 42
SEDR.fix_seed(random_seed)
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

dataset_name = "Gut"
num_clusters = 5
base_path = "../../ST_Datasets"
dataset_path = os.path.join(base_path, dataset_name)

logfile = open("SEDR_Gut_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

print(f"=== SEDR on {dataset_name} (PC=3 and PC=10; k={num_clusters}) ===")
print(f"Device: {device}")
print(f"Dataset path: {dataset_path}")

# --------------------------- Load Data ---------------------------------------
print(f"Loading Visium data from: {dataset_path}")
adata = sc.read_visium(dataset_path, count_file="filtered_feature_bc_matrix.h5", load_images=True)
adata.var_names_make_unique()
print(f"Loaded data: {adata}")

# Subset to selected spots
subset_csv = os.path.join(dataset_path, "gut_df_wt_muscle_rownames.csv")
# If file not found in Gut, try Gut_reduced
if not os.path.exists(subset_csv):
    gut_reduced_path = os.path.join(base_path, "Gut_reduced")
    subset_csv = os.path.join(gut_reduced_path, "gut_df_wt_muscle_rownames.csv")

if os.path.exists(subset_csv):
    spot_ids = pd.read_csv(subset_csv, header=None).squeeze().astype(str)
    adata = adata[adata.obs_names.isin(spot_ids)].copy()
    print(f"✓ Subsetted to {adata.n_obs} spots")
else:
    print(f"Warning: Subsetting file not found at {subset_csv}. Skipping subsetting.")

# --------------------------- Preprocessing ------------------------------------
print("Preprocessing data...")
adata.layers['count'] = adata.X.toarray()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()
print("Preprocessing completed:")
print(adata)

# --------------------------- Construct Graph ----------------------------------
print("Constructing spatial graph...")
graph_dict = SEDR.graph_construction(adata, 12)
print("✓ Spatial graph constructed (k=12).")

# --------------------------- Helper: run one PC setup -----------------------
def run_sedr_for_pcs(n_pcs: int, k_clusters: int = num_clusters):
    assert n_pcs in (3, 10), "This script is set up for 3 or 10 PCs."
    
    print(f"\n--- Running SEDR for n_pcs={n_pcs}, k={k_clusters} ---")
    
    # Compute PCA
    pca_key = f'X_pca_{n_pcs}'
    print(f"Calculating PCA with {n_pcs} components...")
    adata_X = PCA(n_components=n_pcs, random_state=random_seed).fit_transform(adata.layers['count'])
    adata.obsm[pca_key] = adata_X
    print(f"✓ PCA stored in adata.obsm['{pca_key}'] with shape {adata_X.shape}")
    
    # Train SEDR
    sedr_key = f"SEDR_{n_pcs}"
    cluster_key = f"mclust_pca{n_pcs}_{k_clusters}"
    
    print(f"Training SEDR using {pca_key}...")
    sedr_net = SEDR.Sedr(adata.obsm[pca_key], graph_dict, mode='clustering', device=device)
    sedr_net.train_with_dec(N=1)
    sedr_feat, _, _, _ = sedr_net.process()
    adata.obsm[sedr_key] = sedr_feat
    print(f"✓ Stored embedding in adata.obsm['{sedr_key}'] with shape {sedr_feat.shape}")
    
    # Clustering
    SEDR.mclust_R(adata, k_clusters, use_rep=sedr_key, key_added=cluster_key)
    print(f"✓ Clustering done -> obs['{cluster_key}']")
    
    return cluster_key

# --------------------------- Run both PC settings ---------------------------
cluster_keys = []
for pcs in (3, 10):
    key = run_sedr_for_pcs(pcs, k_clusters=num_clusters)
    cluster_keys.append((pcs, key))

# --------------------------- Load and align ground truth ---------------------------
print("\n--- Loading ground truth labels ---")
ari_results = {}
gt_csv_path = os.path.join(dataset_path, "gut_true_labels.csv")
# If file not found in Gut, try Gut_reduced
if not os.path.exists(gt_csv_path):
    gut_reduced_path = os.path.join(base_path, "Gut_reduced")
    gt_csv_path = os.path.join(gut_reduced_path, "gut_true_labels.csv")

if os.path.exists(gt_csv_path):
    try:
        all_true_labels_df = pd.read_csv(gt_csv_path)
        # Swap coordinates as done in STAMarker
        all_true_labels_df.rename(columns={'x': 'y', 'y': 'x'}, inplace=True)
        
        # Align ground truth with experimental data using nearest neighbors
        spatial_coords = pd.DataFrame(adata.obsm['spatial'], columns=['x', 'y'])
        ari_true_labels_df = all_true_labels_df.dropna()
        tree = cKDTree(spatial_coords[['x', 'y']])
        distances, indices = tree.query(ari_true_labels_df[['x', 'y']], k=1)
        
        # Calculate ARI for each PC setting
        for pcs, cluster_key in cluster_keys:
            consensus_aligned = adata.obs[cluster_key].iloc[indices].to_numpy()
            true_aligned = ari_true_labels_df['z'].to_numpy()
            
            gt = pd.Categorical(true_aligned.astype(str)).codes
            pred = pd.Categorical(consensus_aligned.astype(str)).codes
            
            print(f"\nARI calculation for n_pcs={pcs}, k={num_clusters}:")
            print(f"  Number of aligned spots: {len(gt)}")
            print(f"  Unique ground truth values: {sorted(set(true_aligned.astype(str)))}")
            print(f"  Unique predicted clusters: {sorted(set(consensus_aligned.astype(str)))}")
            
            # Create contingency table
            contingency = pd.crosstab(
                pd.Series(true_aligned.astype(str), name="Ground Truth"),
                pd.Series(consensus_aligned.astype(str), name="Predicted")
            )
            print(f"\n  Contingency table (rows=GT, cols=Predicted):")
            print(contingency.to_string())
            
            ari_val = metrics.adjusted_rand_score(gt, pred)
            ari_results[pcs] = ari_val
            print(f"\nARI (n_pcs={pcs}, k={num_clusters}): {ari_val:.6f}")
    except Exception as e:
        print(f"Error loading or aligning ground truth: {e}")
        import traceback
        traceback.print_exc()
else:
    print(f"Ground truth file not found at: {gt_csv_path}")

# --------------------------- Save combined results --------------------------
print("\n--- Saving combined results ---")
combined = pd.DataFrame({
    "barcode": adata.obs_names.astype(str).values,
    "x": adata.obsm['spatial'][:, 0],
    "y": adata.obsm['spatial'][:, 1],
})

for pcs, key in cluster_keys:
    combined[f"label_pca{pcs}_k{num_clusters}"] = adata.obs[key].astype(str).values

# Add ground truth if available
if os.path.exists(gt_csv_path):
    try:
        # Merge ground truth by coordinates
        combined['x_rounded'] = combined['x'].round().astype(int)
        combined['y_rounded'] = combined['y'].round().astype(int)
        all_true_labels_df['x'] = all_true_labels_df['x'].round().astype(int)
        all_true_labels_df['y'] = all_true_labels_df['y'].round().astype(int)
        
        merged = pd.merge(
            combined,
            all_true_labels_df[['x', 'y', 'z']],
            left_on=['x_rounded', 'y_rounded'],
            right_on=['x', 'y'],
            how='left'
        )
        combined['true_label'] = merged['z']
        combined = combined.drop(columns=['x_rounded', 'y_rounded'])
    except Exception as e:
        print(f"Warning: Could not merge ground truth: {e}")

combined_out = "SEDR_Gut.csv"
combined.to_csv(combined_out, index=False)
print(f"✓ Saved combined results: {combined_out}")

# Save ARI summary if computed
if ari_results:
    print("\n--- ARI Summary ---")
    ari_lines = []
    for pcs, ari_val in ari_results.items():
        line = f"pca={pcs}, k={num_clusters} -> ARI={ari_val:.6f}"
        print(line)
        ari_lines.append(line)
    with open("gut_sedr_ari.txt", "w") as f:
        for line in ari_lines:
            f.write(line + "\n")
    print("✓ Saved ARI summary to gut_sedr_ari.txt")

print(f"\n=== Done: SEDR on {dataset_name} for PC=3 and PC=10 with k={num_clusters} ===")
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__
logfile.close()

