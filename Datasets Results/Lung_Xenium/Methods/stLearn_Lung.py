#!/usr/bin/env python3
"""
stLearn clustering analysis for Lung dataset with histology image processing.
"""
import numpy as np
import pandas as pd
from pathlib import Path
import scanpy as sc
import stlearn as st
from sklearn import metrics
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import warnings
import random
import os
import gc

warnings.filterwarnings("ignore")
st.settings.set_figure_params(dpi=120)

# Set seeds for reproducibility
seed = 42
np.random.seed(seed)
random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

# Dataset configuration
DATASET_NAME = "Lung"
K = 6  # True number of clusters
H5AD_PATH = "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/Data_VUILD96MF/Lung_xenium_processed.h5ad"
OUT_DIR = Path(".")
TILE_DIR = Path("/tmp/tiles_lung")
TILE_DIR.mkdir(parents=True, exist_ok=True)

PCS_LIST = [3, 10]
N_COMPS = max(PCS_LIST)

print(f"=== Processing {DATASET_NAME} ===")
print(f"H5AD path: {H5AD_PATH}")

# Load data
print("\nLoading data...")
adata = sc.read_h5ad(H5AD_PATH)
adata.var_names_make_unique()

# Ensure spatial coordinates are in obsm
if "spatial" not in adata.obsm:
    if "x" in adata.obs.columns and "y" in adata.obs.columns:
        adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()
    elif "X" in adata.obs.columns and "Y" in adata.obs.columns:
        adata.obsm["spatial"] = adata.obs[["X", "Y"]].to_numpy()
    elif "x_centroid" in adata.obs.columns and "y_centroid" in adata.obs.columns:
        adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
    else:
        print("Available columns:", adata.obs.columns.tolist())
        print("Available obsm keys:", list(adata.obsm.keys()))
        raise ValueError("Could not find spatial coordinates")

# Store original counts (only if needed, to save memory)
# adata.layers["counts"] = adata.X.copy()  # Commented out to save memory

# Create dummy Visium-like fields (required by stLearn)
xy = adata.obsm["spatial"]
adata.obs["imagerow"] = np.round(xy[:, 1]).astype(int)
adata.obs["imagecol"] = np.round(xy[:, 0]).astype(int)
adata.obs["array_row"] = adata.obs["imagerow"]
adata.obs["array_col"] = adata.obs["imagecol"]

print(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
# Calculate memory usage (handle both sparse and dense matrices)
if hasattr(adata.X, 'data'):
    # Sparse matrix
    mem_gb = (adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes) / 1024**3
else:
    # Dense matrix
    mem_gb = adata.X.nbytes / 1024**3
print(f"Data shape: {adata.X.shape}, Memory usage: {mem_gb:.2f} GB")

# Pre-processing for gene count table
print("\nPre-processing gene count table...")
st.pp.filter_genes(adata, min_cells=1)
st.pp.normalize_total(adata)
st.pp.log1p(adata)

# Pre-processing for spot image (histology)
print("\nPre-processing histology images...")
try:
    st.pp.tiling(adata, out_path=str(TILE_DIR))
    st.pp.extract_feature(adata)
    print("Histology features extracted successfully")
except Exception as e:
    print(f"Warning: Could not extract histology features: {e}")
    print("Creating dummy morphology features...")
    if "X_morphology" not in adata.obsm:
        adata.obsm["X_morphology"] = np.zeros((adata.n_obs, 8), dtype=float)

# Run initial PCA (needed by SME)
print(f"\nRunning initial PCA with {N_COMPS} components...")
st.em.run_pca(adata, n_comps=N_COMPS)

# Apply SME normalization
print("\nApplying stSME normalization...")
adata_sme = adata.copy()
# Free up memory from original adata after copy
del adata
gc.collect()

# Apply SME normalize on the log-transformed data
st.spatial.SME.SME_normalize(adata_sme, use_data="raw", platform="Visium", weights="weights_matrix_pd_gd")

# Guard against NaNs/Infs from SME
adata_sme.obsm["raw_SME_normalized"] = np.nan_to_num(
    adata_sme.obsm["raw_SME_normalized"],
    nan=0.0,
    posinf=0.0,
    neginf=0.0
)

# Set SME-normalized data as main matrix
adata_sme.X = adata_sme.obsm['raw_SME_normalized']

# Scale and run PCA on SME-normalized data
print("\nScaling SME-normalized data...")
st.pp.scale(adata_sme)
adata_sme.X = np.nan_to_num(adata_sme.X, nan=0.0, posinf=0.0, neginf=0.0)

print(f"Running final PCA with {N_COMPS} components on SME-normalized data...")
st.em.run_pca(adata_sme, n_comps=N_COMPS)

# Run clustering for different PC counts
print("\nRunning k-means clustering...")
for n_pcs in PCS_LIST:
    print(f"  Clustering with {n_pcs} PCs, k={K}...")
    adata_sme.obsm["X_pca_subset"] = adata_sme.obsm["X_pca"][:, :n_pcs].copy()
    st.tl.clustering.kmeans(
        adata_sme,
        n_clusters=K,
        use_data="X_pca_subset",
        key_added=f"label_{n_pcs}PCs",
        algorithm="lloyd"
    )
    print(f"    Clustering complete for {n_pcs} PCs")

# Save combined results
print("\nSaving results...")
xy = adata_sme.obsm["spatial"]
results_df = pd.DataFrame({
    "ID": adata_sme.obs_names,
    "x": xy[:, 0],
    "y": xy[:, 1],
    f"label_3PCs_k{K}": adata_sme.obs["label_3PCs"],
    f"label_10PCs_k{K}": adata_sme.obs["label_10PCs"]
})

csv_path = OUT_DIR / f"stLearn_{DATASET_NAME}.csv"
results_df.to_csv(csv_path, index=False)
print(f"Saved results to: {csv_path}")

# Find ground truth column
ground_truth_col = None
for gt_col in ['Annotation_Type', 'ground_truth', 'true_labels', 'label', 'cluster', 'cell_type']:
    if gt_col in adata_sme.obs.columns:
        ground_truth_col = gt_col
        break

# Load ground truth from CSV if not found in obs
loaded_from_csv = False
if not ground_truth_col:
    # Try multiple possible CSV locations
    h5ad_dir = Path(H5AD_PATH).parent
    possible_csv_paths = [
        h5ad_dir / "lung_true_labels.csv",
        h5ad_dir / "Lung_true_labels.csv",
        h5ad_dir / "true_labels.csv",
        h5ad_dir.parent / "lung_true_labels.csv",
        h5ad_dir.parent / "Lung_true_labels.csv",
    ]
    
    csv_path = None
    for path in possible_csv_paths:
        if path.exists():
            csv_path = path
            break
    
    if csv_path:
        print(f"\nLoading ground truth from CSV: {csv_path}")
        all_true_labels_df = pd.read_csv(csv_path)
        
        # Handle different column name conventions
        # Check if we need to swap x/y (some datasets have them swapped)
        if 'x' in all_true_labels_df.columns and 'y' in all_true_labels_df.columns:
            # Check which column contains labels (could be 'z', 'label', 'cluster', 'cell_type', etc.)
            label_col = None
            for col in ['z', 'label', 'cluster', 'cell_type', 'type', 'ground_truth']:
                if col in all_true_labels_df.columns:
                    label_col = col
                    break
            
            if label_col:
                spatial_coords = pd.DataFrame(adata_sme.obsm['spatial'], columns=['x', 'y'])
                ari_true_labels_df = all_true_labels_df.dropna()
                
                if len(ari_true_labels_df) > 0:
                    tree = cKDTree(spatial_coords[['x', 'y']])
                    distances, indices = tree.query(ari_true_labels_df[['x', 'y']], k=1)
                    
                    adata_sme.obs['true_labels'] = pd.NA
                    aligned_labels = ari_true_labels_df[label_col].astype(str).values
                    adata_sme.obs.iloc[indices, adata_sme.obs.columns.get_loc('true_labels')] = aligned_labels
                    
                    ground_truth_col = 'true_labels'
                    loaded_from_csv = True
                    print(f"Aligned {len(ari_true_labels_df)} ground truth labels from column '{label_col}'")
            else:
                print("Warning: Found CSV but could not identify label column")
        else:
            print(f"Warning: CSV file found but missing 'x' and 'y' columns. Available columns: {all_true_labels_df.columns.tolist()}")

# Compute ARI scores
print("\nComputing ARI scores...")
ari_results = {}
if ground_truth_col:
    print(f"Using '{ground_truth_col}' as ground truth")
    if loaded_from_csv and ground_truth_col == 'true_labels':
        aligned_mask = ~adata_sme.obs[ground_truth_col].isna()
        if aligned_mask.sum() > 0:
            gt = pd.Categorical(adata_sme.obs.loc[aligned_mask, ground_truth_col].astype(str)).codes
            for n_pcs in PCS_LIST:
                pr = pd.Categorical(adata_sme.obs.loc[aligned_mask, f"label_{n_pcs}PCs"]).codes
                ari = metrics.adjusted_rand_score(gt, pr)
                ari_results[n_pcs] = ari
                print(f"  ARI ({n_pcs} PCs, k={K}): {ari:.6f}")
            if aligned_mask.sum() < len(adata_sme):
                print(f"Note: Calculated ARI on {aligned_mask.sum()} out of {len(adata_sme)} observations")
        else:
            print("Warning: No aligned ground truth labels found (all are NA)")
    else:
        valid_mask = ~adata_sme.obs[ground_truth_col].isna()
        if valid_mask.sum() > 0:
            gt = pd.Categorical(adata_sme.obs.loc[valid_mask, ground_truth_col].astype(str)).codes
            for n_pcs in PCS_LIST:
                pr = pd.Categorical(adata_sme.obs.loc[valid_mask, f"label_{n_pcs}PCs"]).codes
                ari = metrics.adjusted_rand_score(gt, pr)
                ari_results[n_pcs] = ari
                print(f"  ARI ({n_pcs} PCs, k={K}): {ari:.6f}")
            if valid_mask.sum() < len(adata_sme):
                print(f"Note: Calculated ARI on {valid_mask.sum()} out of {len(adata_sme)} observations")
        else:
            print("Warning: No valid ground truth labels found (all are NA)")
else:
    print("Warning: No ground truth column found")
    print(f"Available obs columns: {list(adata_sme.obs.columns)}")
    print("Please check if ground truth CSV file exists or if ground truth is stored under a different column name")

# Save ARI summary
ari_path = OUT_DIR / f"stLearn_{DATASET_NAME}_ARI.txt"
with open(ari_path, "w") as f:
    f.write(f"stSME-based clustering results for {DATASET_NAME}\n")
    f.write("=" * 50 + "\n")
    if ground_truth_col:
        f.write(f"Ground truth column: {ground_truth_col}\n")
    f.write(f"Number of clusters: {K}\n")
    f.write("=" * 50 + "\n")
    for n_pcs in PCS_LIST:
        if n_pcs in ari_results:
            f.write(f"{n_pcs}PCs: ARI = {ari_results[n_pcs]:.6f}\n")
        else:
            f.write(f"{n_pcs}PCs: ARI calculation skipped\n")
print(f"Saved ARI summary to: {ari_path}")

# Create plot with 3 columns (3PCs, 10PCs, ground truth)
if "spatial" in adata_sme.obsm:
    print("\nCreating visualization plot...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    spatial_coords = adata_sme.obsm['spatial']
    
    # Plot 3PCs clustering
    scatter1 = axes[0].scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                               c=pd.Categorical(adata_sme.obs["label_3PCs"]).codes,
                               cmap='tab20', s=5)
    axes[0].set_title(f'{DATASET_NAME} - 3 PCs, k={K}')
    axes[0].set_xlabel('X coordinate')
    axes[0].set_ylabel('Y coordinate')
    axes[0].set_aspect('equal', adjustable='box')
    
    # Plot 10PCs clustering
    scatter2 = axes[1].scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                               c=pd.Categorical(adata_sme.obs["label_10PCs"]).codes,
                               cmap='tab20', s=5)
    axes[1].set_title(f'{DATASET_NAME} - 10 PCs, k={K}')
    axes[1].set_xlabel('X coordinate')
    axes[1].set_ylabel('Y coordinate')
    axes[1].set_aspect('equal', adjustable='box')
    
    # Plot ground truth
    if ground_truth_col:
        scatter3 = axes[2].scatter(spatial_coords[:, 0], spatial_coords[:, 1],
                                   c=pd.Categorical(adata_sme.obs[ground_truth_col]).codes,
                                   cmap='tab20', s=5)
        axes[2].set_title(f'{DATASET_NAME} - Ground Truth')
        axes[2].set_xlabel('X coordinate')
        axes[2].set_ylabel('Y coordinate')
        axes[2].set_aspect('equal', adjustable='box')
    else:
        axes[2].text(0.5, 0.5, 'No ground truth\navailable',
                    ha='center', va='center', transform=axes[2].transAxes)
        axes[2].set_title(f'{DATASET_NAME} - Ground Truth (N/A)')
    
    plt.tight_layout()
    plot_path = OUT_DIR / f"stLearn_{DATASET_NAME}_plot.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to: {plot_path}")

print(f"\n✓ {DATASET_NAME} analysis complete!")

