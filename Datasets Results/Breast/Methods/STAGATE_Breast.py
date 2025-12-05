#!/usr/bin/env python3
import sys
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
sys.path.insert(0, os.path.join(project_root, 'STAGATE_local_code'))

os.environ['CUDA_HOME'] = '/sw/eb/sw/CUDA/11.7.0'
os.environ['PATH'] = f"{os.environ['CUDA_HOME']}/bin:{os.environ['PATH']}"
os.environ['LD_LIBRARY_PATH'] = f"{os.environ['CUDA_HOME']}/lib64:{os.environ.get('LD_LIBRARY_PATH', '')}"

os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a'
os.environ['R_LIBS'] = f"/scratch/user/varogovchenko/Rlibs:/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R/library"
os.environ['PATH'] = f"{os.environ['R_HOME']}/bin:{os.environ['PATH']}"

logfile = open("STAGATE_Breast_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import torch
import tensorflow as tf
from sklearn import metrics
from scipy.spatial import cKDTree
import random

# Set seeds for reproducibility
seed = 42
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
tf.random.set_seed(seed)
if torch.cuda.is_available():
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

tf.compat.v1.disable_eager_execution()
physical_devices = tf.config.list_physical_devices('GPU')
if physical_devices:
    try:
        for device in physical_devices:
            tf.config.experimental.set_memory_growth(device, True)
    except RuntimeError as e:
        print(f"GPU config error: {e}")

from STAGATE.utils import Cal_Spatial_Net, Stats_Spatial_Net, mclust_R
from STAGATE.Train_STAGATE import train_STAGATE

dataset_name = "Breast"
num_clusters = 5
base_path = "../../ST_Datasets"
dataset_path = os.path.join(base_path, dataset_name)

print(f"Processing {dataset_name}...")
print(f"Path: {dataset_path}")

adata = sc.read_visium(dataset_path, count_file="filtered_feature_bc_matrix.h5", load_images=True)
adata.var_names_make_unique()
print(adata)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()
print(adata)

print("\nRunning STAGATE...")
Cal_Spatial_Net(adata, rad_cutoff=150)
Stats_Spatial_Net(adata)
adata = train_STAGATE(adata, alpha=0, hidden_dims=[512, 30])
sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)
adata = mclust_R(adata, used_obsm='STAGATE', num_cluster=num_clusters)
adata.obs['mclust'] = adata.obs['mclust'].astype(str)

df = pd.DataFrame({
    'barcode': adata.obs_names,
    'mclust': adata.obs['mclust'].values
})
df[['x', 'y']] = adata.obsm['spatial']

base_dir = os.path.dirname(os.path.abspath(__file__))
results_file = os.path.join(base_dir, "STAGATE_Breast.csv")
df.to_csv(results_file, index=False)
print(f"Saved results to: {results_file}")

ari = None
ground_truth_col = None

for gt_col in ['ground_truth', 'true_labels', 'label', 'cluster', 'cell_type']:
    if gt_col in adata.obs.columns:
        ground_truth_col = gt_col
        break

true_labels_csv = os.path.join(dataset_path, "breast_true_labels.csv")
loaded_from_csv = False

if not ground_truth_col:
    if os.path.exists(true_labels_csv):
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

if ground_truth_col:
    print(f"\nComputing ARI (using '{ground_truth_col}')...")
    
    if loaded_from_csv and ground_truth_col == 'true_labels':
        aligned_mask = ~adata.obs[ground_truth_col].isna()
        if aligned_mask.sum() > 0:
            gt = pd.Categorical(adata.obs.loc[aligned_mask, ground_truth_col].astype(str)).codes
            pr = pd.Categorical(adata.obs.loc[aligned_mask, 'mclust']).codes
            
            ari = metrics.adjusted_rand_score(gt, pr)
            print(f"ARI: {ari:.6f}")
    else:
        valid_mask = ~adata.obs[ground_truth_col].isna()
        if valid_mask.sum() > 0:
            gt = pd.Categorical(adata.obs.loc[valid_mask, ground_truth_col].astype(str)).codes
            pr = pd.Categorical(adata.obs.loc[valid_mask, 'mclust']).codes
            
            ari = metrics.adjusted_rand_score(gt, pr)
            print(f"ARI: {ari:.6f}")
            if valid_mask.sum() < len(adata):
                print(f"Note: ARI calculated on {valid_mask.sum()}/{len(adata)} observations (excluded {len(adata) - valid_mask.sum()} NA)")
    
    ari_file = os.path.join(base_dir, "STAGATE_Breast_ARI.txt")
    with open(ari_file, "w") as f:
        f.write(f"STAGATE clustering results for {dataset_name}\n")
        f.write("="*50 + "\n")
        f.write(f"Ground truth column: {ground_truth_col}\n")
        f.write(f"Number of clusters: {num_clusters}\n")
        f.write("="*50 + "\n")
        if ari is not None:
            f.write(f"ARI: {ari:.6f}\n")
        else:
            f.write("ARI calculation skipped: No valid ground truth labels found.\n")
    print(f"Saved ARI to: {ari_file}")
else:
    print("No ground truth column found - skipping ARI")
    ari_file = os.path.join(base_dir, "STAGATE_Breast_ARI.txt")
    with open(ari_file, "w") as f:
        f.write(f"STAGATE clustering results for {dataset_name}\n")
        f.write("="*50 + "\n")
        f.write("No ground truth column found — ARI computation skipped.\n")

fig, axs = plt.subplots(1, 2, figsize=(12, 5))

sc.pl.spatial(adata, color='mclust', ax=axs[0], show=False, 
              title=f'{dataset_name} - STAGATE', size=1.5)

if ground_truth_col and ground_truth_col in adata.obs.columns:
    if ground_truth_col == 'true_labels':
        aligned_mask = ~adata.obs['true_labels'].isna()
        if aligned_mask.sum() > 0:
            temp_adata = adata[aligned_mask].copy()
            sc.pl.spatial(temp_adata, color=ground_truth_col, ax=axs[1], 
                         show=False, title=f'{dataset_name} - True Labels', size=1.5)
        else:
            axs[1].text(0.5, 0.5, 'No aligned ground truth\nlabels available', 
                       ha='center', va='center', transform=axs[1].transAxes)
            axs[1].set_title(f'{dataset_name} - True Labels (N/A)')
    else:
        sc.pl.spatial(adata, color=ground_truth_col, ax=axs[1], show=False, 
                     title=f'{dataset_name} - True Labels', size=1.5)
else:
    axs[1].text(0.5, 0.5, 'No ground truth\navailable', 
               ha='center', va='center', transform=axs[1].transAxes)
    axs[1].set_title(f'{dataset_name} - True Labels (N/A)')

plt.tight_layout()
plot_file = os.path.join(base_dir, "STAGATE_Breast_plot.png")
plt.savefig(plot_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved plot to: {plot_file}")

print(f"\nFinished {dataset_name}")
logfile.close()
