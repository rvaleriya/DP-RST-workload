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
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import cv2

# =========================
# Logging
# =========================
LOGFILE = "SpaGCN_Breast_histology_output.log"
logfile = open(LOGFILE, "w")
sys.stdout = logfile
sys.stderr = logfile
print(f"--- Script Start: Logging to {LOGFILE} ---")

# =========================
# Config
# =========================
BASE_DATA_DIR = Path("/scratch/user/varogovchenko/ST_Datasets")
DATASET_NAME = "Breast"
DATASET_PATH = BASE_DATA_DIR / DATASET_NAME
COUNT_FILE = "filtered_feature_bc_matrix.h5"
IMAGE_FILE = "spatial/tissue_hires_image.png"
OUTPUT_CSV = "SpaGCN_Breast_histology.csv"
OUTPUT_ARI_TXT = "SpaGCN_Breast_histology_ARI.txt"
OUTPUT_PLOT = "SpaGCN_Breast_histology_plot.png"

N_CLUSTERS = 5
NUM_PCS_LIST = [3, 10]
P_PARAM = 0.5
SEED = 42
ALPHA = 1
BETA = 49

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
# Load Visium data and histology image
# =========================
print(f"Loading Visium data from: {DATASET_PATH}")
adata = sc.read_visium(str(DATASET_PATH), count_file=COUNT_FILE, load_images=True)
adata.var_names_make_unique()
print(adata)

# Load histology image
image_path = DATASET_PATH / IMAGE_FILE
if not image_path.exists():
    image_path = DATASET_PATH / "spatial/tissue_tiff_image.tif"
print(f"Loading histology image from: {image_path}")
img = cv2.imread(str(image_path))
if img is None:
    print(f"Warning: Could not load image from {image_path}. Proceeding without histology image.")
    img = None
else:
    print(f"Loaded histology image with shape: {img.shape}")

# Coordinates
coords = adata.obsm["spatial"]
x_array = adata.obs["array_col"].tolist()
y_array = adata.obs["array_row"].tolist()

# Scale coordinates to match image dimensions
# Visium stores coordinates in full-resolution space, but tissue_hires_image.png is scaled
if img is not None and 'spatial' in adata.uns:
    # Get scaling factors from Visium metadata
    sample_id = list(adata.uns['spatial'].keys())[0]
    spatial_info = adata.uns['spatial'][sample_id]
    
    # Get full resolution and scaled image dimensions
    if 'scalefactors' in spatial_info:
        scalefactors = spatial_info['scalefactors']
        # tissue_hires_image.png scale factor
        scale_factor = scalefactors.get('tissue_hires_scalef', None)
        
        if scale_factor is not None:
            print(f"Using Visium scale factor: {scale_factor}")
            x_pixel_full = coords[:, 0]
            y_pixel_full = coords[:, 1]
            # Convert to integers for image indexing
            x_pixel = (x_pixel_full * scale_factor).astype(int).tolist()
            y_pixel = (y_pixel_full * scale_factor).astype(int).tolist()
        else:
            # Fallback: calculate scale factor from image dimensions
            img_height, img_width = img.shape[0], img.shape[1]
            x_range = coords[:, 0].max() - coords[:, 0].min()
            y_range = coords[:, 1].max() - coords[:, 1].min()
            scale_x = img_width / x_range if x_range > 0 else 1.0
            scale_y = img_height / y_range if y_range > 0 else 1.0
            scale_factor = min(scale_x, scale_y) * 0.9  # Use 90% to ensure within bounds
            print(f"Calculated scale factor: {scale_factor} (image: {img_width}x{img_height}, coords: {x_range:.1f}x{y_range:.1f})")
            x_pixel_full = coords[:, 0]
            y_pixel_full = coords[:, 1]
            # Convert to integers for image indexing
            x_pixel = ((x_pixel_full - coords[:, 0].min()) * scale_factor).astype(int).tolist()
            y_pixel = ((y_pixel_full - coords[:, 1].min()) * scale_factor).astype(int).tolist()
    else:
        # No scalefactors found, try to infer from image and coordinate ranges
        img_height, img_width = img.shape[0], img.shape[1]
        x_range = coords[:, 0].max() - coords[:, 0].min()
        y_range = coords[:, 1].max() - coords[:, 1].min()
        scale_x = img_width / x_range if x_range > 0 else 1.0
        scale_y = img_height / y_range if y_range > 0 else 1.0
        scale_factor = min(scale_x, scale_y) * 0.9
        print(f"Inferred scale factor: {scale_factor} (image: {img_width}x{img_height}, coords: {x_range:.1f}x{y_range:.1f})")
        x_pixel_full = coords[:, 0]
        y_pixel_full = coords[:, 1]
        # Convert to integers for image indexing
        x_pixel = ((x_pixel_full - coords[:, 0].min()) * scale_factor).astype(int).tolist()
        y_pixel = ((y_pixel_full - coords[:, 1].min()) * scale_factor).astype(int).tolist()
else:
    # Convert to integers for image indexing
    x_pixel = coords[:, 0].astype(int).tolist()
    y_pixel = coords[:, 1].astype(int).tolist()

adata.obs["x_array"] = x_array
adata.obs["y_array"] = y_array
adata.obs["x_pixel"] = x_pixel
adata.obs["y_pixel"] = y_pixel

# Check coordinate ranges for debugging
if img is not None:
    print(f"Image shape: {img.shape} (height, width, channels)")
    print(f"Scaled coordinate ranges: x_pixel [{min(x_pixel):.1f}, {max(x_pixel):.1f}], y_pixel [{min(y_pixel):.1f}, {max(y_pixel):.1f}]")

# =========================
# Preprocessing
# =========================
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
max_pcs = max(NUM_PCS_LIST)
sc.tl.pca(adata, n_comps=max_pcs)
print(f"✓ Computed {max_pcs} PCs")

# =========================
# Adjacency with histology
# =========================
use_histology_adj = False
if img is not None:
    print(f"Calculating adjacency with histology image (alpha={ALPHA}, beta={BETA})...")
    # Try extracting color first to test if coordinates work
    try:
        # Test color extraction for a few spots
        beta_half = round(BETA / 2)
        test_indices = [0, len(x_pixel)//2, len(x_pixel)-1] if len(x_pixel) > 2 else [0]
        color_test_passed = True
        for idx in test_indices:
            max_x, max_y = img.shape[0], img.shape[1]
            x_coord, y_coord = int(x_pixel[idx]), int(y_pixel[idx])
            if x_coord < 0 or x_coord >= max_x or y_coord < 0 or y_coord >= max_y:
                print(f"Warning: Coordinates out of bounds at index {idx}: ({x_coord}, {y_coord}) vs image size ({max_x}, {max_y})")
                color_test_passed = False
                break
            nbs = img[max(0, x_coord-beta_half):min(max_x, x_coord+beta_half+1),
                      max(0, y_coord-beta_half):min(max_y, y_coord+beta_half+1)]
            if nbs.size == 0:
                print(f"Warning: Empty neighborhood at index {idx}")
                color_test_passed = False
                break
        
        if color_test_passed:
            adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, 
                                           image=img, beta=BETA, alpha=ALPHA, histology=True)
            # Check if adjacency matrix contains NaN values
            if np.isnan(adj).any() or np.isinf(adj).any():
                raise RuntimeError("Histology-based adjacency matrix contains NaN/Inf values. Check coordinate scaling and image alignment.")
            use_histology_adj = True
            print("✓ Successfully calculated histology-based adjacency matrix")
        else:
            raise RuntimeError("Color extraction test failed. Coordinates are out of image bounds. Check coordinate scaling.")
    except Exception as e:
        raise RuntimeError(f"Failed to calculate histology-based adjacency: {e}")
else:
    raise RuntimeError("No histology image available. Cannot proceed without histology image.")

# =========================
# Search l for p
# =========================
print(f"Searching l for p={P_PARAM} ...")
l = spg.search_l(P_PARAM, adj)
if l is None:
    print("Warning: search_l returned None. Falling back to spatial-only adjacency for l calculation.")
    adj_spatial = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    l = spg.search_l(P_PARAM, adj_spatial)
    if l is None:
        raise RuntimeError("Failed to calculate l parameter even with spatial-only adjacency. Check adjacency matrix.")
    print(f"Using spatial-only l = {l}")
    # If histology adj had issues, use spatial adj
    if not use_histology_adj:
        adj = adj_spatial
else:
    print(f"Recommended l = {l}")

# =========================
# Load ground truth labels
# =========================
ground_truth_col = None
for gt_col in ['ground_truth', 'true_labels', 'label', 'cluster', 'cell_type']:
    if gt_col in adata.obs.columns:
        ground_truth_col = gt_col
        break

true_labels_csv = DATASET_PATH / "breast_true_labels.csv"
loaded_from_csv = False

if not ground_truth_col:
    if true_labels_csv.exists():
        print(f"Loading ground truth from CSV: {true_labels_csv}")
        all_true_labels_df = pd.read_csv(true_labels_csv)
        all_true_labels_df.rename(columns={'x': 'y', 'y': 'x'}, inplace=True)
        
        spatial_coords = pd.DataFrame(adata.obsm['spatial'], columns=['x', 'y'])
        ari_true_labels_df = all_true_labels_df.dropna()
        
        if len(ari_true_labels_df) > 0:
            tree = cKDTree(spatial_coords[['x', 'y']])
            distances, indices = tree.query(ari_true_labels_df[['x', 'y']], k=1)
            
            adata.obs['true_labels'] = pd.NA
            aligned_labels = ari_true_labels_df['z'].astype(str).values
            adata.obs.iloc[indices, adata.obs.columns.get_loc('true_labels')] = aligned_labels
            
            ground_truth_col = 'true_labels'
            loaded_from_csv = True
            print(f"Aligned {len(ari_true_labels_df)} ground truth labels")

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
    print(f"Running SpaGCN with {num_pcs} PCs (k={N_CLUSTERS}) - WITH HISTOLOGY")
    print("="*60)

    adata_run = adata.copy()
    adata_run.obsm["X_pca"] = adata_run.obsm["X_pca"][:, :num_pcs]

    r_seed = t_seed = n_seed = SEED
    print(f"Searching resolution for n_clusters={N_CLUSTERS} ...")
    res = spg.search_res(
        adata_run, adj, l, N_CLUSTERS,
        start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20,
        r_seed=r_seed, t_seed=t_seed, n_seed=n_seed, num_pcs=num_pcs
    )
    print("Recommended res =", res)

    clf = spg.SpaGCN()
    if l is None:
        raise RuntimeError(f"l parameter is None for {num_pcs} PCs. Cannot proceed with training.")
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

    print("Refining clusters (shape='hexagon')...")
    # Use spatial-only adjacency for refinement (as per tutorial)
    adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    refined = spg.refine(
        sample_id=adata_run.obs.index.tolist(),
        pred=adata_run.obs["pred"].tolist(),
        dis=adj_2d,
        shape="hexagon"
    )
    adata_run.obs["refined_pred"] = pd.Categorical(refined)
    print("Refined prediction stored in adata_run.obs['refined_pred']")

    col_init = f"pred_{num_pcs}PCs_k{N_CLUSTERS}"
    col_ref  = f"refined_pred_{num_pcs}PCs_k{N_CLUSTERS}"
    results[col_init] = adata_run.obs["pred"].astype(str).values
    results[col_ref]  = adata_run.obs["refined_pred"].astype(str).values

    # ARI calculation
    if ground_truth_col:
        if loaded_from_csv and ground_truth_col == 'true_labels':
            aligned_mask = ~adata_run.obs[ground_truth_col].isna()
            if aligned_mask.sum() > 0:
                gt = adata_run.obs.loc[aligned_mask, ground_truth_col].astype(str).tolist()
                pr_i = adata_run.obs.loc[aligned_mask, "pred"].astype(str).tolist()
                pr_r = adata_run.obs.loc[aligned_mask, "refined_pred"].astype(str).tolist()
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
        else:
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
        f.write("SpaGCN on Breast (k=5) WITH HISTOLOGY — ARI summary\n")
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
print("Reloading data with images for plotting...")
# Always reload with images to ensure image is available
adata_plot_base = sc.read_visium(str(DATASET_PATH), count_file=COUNT_FILE, load_images=True)
adata_plot_base.var_names_make_unique()

# Transfer clustering results and ground truth labels
results_indexed = results.set_index("barcode")
adata_plot_base.obs["refined_3PCs"] = results_indexed.loc[adata_plot_base.obs_names, "refined_pred_3PCs_k5"].astype('category')
adata_plot_base.obs["refined_10PCs"] = results_indexed.loc[adata_plot_base.obs_names, "refined_pred_10PCs_k5"].astype('category')
if ground_truth_col and ground_truth_col in adata.obs.columns:
    adata_plot_base.obs[ground_truth_col] = adata.obs.loc[adata_plot_base.obs_names, ground_truth_col]

# Check if image is available
has_image = 'spatial' in adata_plot_base.uns and isinstance(adata_plot_base.uns['spatial'], dict)
if has_image:
    sample_id = list(adata_plot_base.uns['spatial'].keys())[0]
    has_hires = 'images' in adata_plot_base.uns['spatial'][sample_id] and 'hires' in adata_plot_base.uns['spatial'][sample_id]['images']
    print(f"Image check: spatial dict exists={has_image}, sample_id={sample_id}, hires image={has_hires}")
else:
    print("Warning: No spatial image found in adata.uns['spatial']")

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 3PCs results
sc.pl.spatial(adata_plot_base, color="refined_3PCs", ax=axes[0], show=False, 
              title=f'Breast - 3PCs (k={N_CLUSTERS}) [Histology]', spot_size=120)

# 10PCs results
sc.pl.spatial(adata_plot_base, color="refined_10PCs", ax=axes[1], show=False, 
              title=f'Breast - 10PCs (k={N_CLUSTERS}) [Histology]', spot_size=120)

# True labels
if ground_truth_col and ground_truth_col in adata_plot_base.obs.columns:
    if ground_truth_col == 'true_labels':
        aligned_mask = ~adata_plot_base.obs['true_labels'].isna()
        if aligned_mask.sum() > 0:
            temp_adata = adata_plot_base[aligned_mask].copy()
            sc.pl.spatial(temp_adata, color=ground_truth_col, ax=axes[2], 
                         show=False, title='Breast - True Labels', spot_size=120)
        else:
            axes[2].text(0.5, 0.5, 'No aligned ground truth\nlabels available', 
                       ha='center', va='center', transform=axes[2].transAxes)
            axes[2].set_title('Breast - True Labels (N/A)')
    else:
        sc.pl.spatial(adata_plot_base, color=ground_truth_col, ax=axes[2], show=False, 
                     title='Breast - True Labels', spot_size=120)
else:
    axes[2].text(0.5, 0.5, 'No ground truth\navailable', 
               ha='center', va='center', transform=axes[2].transAxes)
    axes[2].set_title('Breast - True Labels (N/A)')

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved plot to {OUTPUT_PLOT}")

print("\nDone.")
logfile.close()

