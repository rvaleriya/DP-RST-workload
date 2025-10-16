#!/usr/bin/env python3
import os
import sys
import random
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'  # disable GPU if needed
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'   # suppress TF messages

# CUDA paths (safe even if GPU disabled)
os.environ['CUDA_HOME'] = '/sw/eb/sw/CUDA/11.7.0'
os.environ['PATH'] = f"{os.environ['CUDA_HOME']}/bin:{os.environ['PATH']}"
os.environ['LD_LIBRARY_PATH'] = f"{os.environ['CUDA_HOME']}/lib64:{os.environ.get('LD_LIBRARY_PATH', '')}"

# R paths for rpy2 + mclust
os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a'
os.environ['R_LIBS'] = f"/scratch/user/varogovchenko/Rlibs:/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R/library"
os.environ['PATH'] = f"{os.environ['R_HOME']}/bin:{os.environ['PATH']}"

# Logging
logfile = open("STAGATE_STARmap_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import pandas as pd
import torch
import tensorflow as tf
from sklearn import metrics

# STAGATE imports
from STAGATE.utils import Cal_Spatial_Net, Stats_Spatial_Net, mclust_R
from STAGATE.Train_STAGATE import train_STAGATE

# ----------------- Config -----------------
H5AD = "/scratch/user/varogovchenko/BASTION_HPRC/STARmap/STARmap_20180505_BY3_1k_20251008011714.h5ad"
OUTPUT_CSV = "STAGATE_STARmap_jointPCA.csv"
OUTPUT_ARI_TXT = "STAGATE_STARmap_ARI.txt"
N_CLUSTERS = 7
K = N_CLUSTERS

# TensorFlow configuration
tf.compat.v1.disable_eager_execution()
physical_devices = tf.config.list_physical_devices('GPU')
if physical_devices:
    try:
        for device in physical_devices:
            tf.config.experimental.set_memory_growth(device, True)
        print(f"Found {len(physical_devices)} GPU(s)")
    except RuntimeError as e:
        print(f"Error configuring GPU: {e}")
else:
    print("No GPU devices found")

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

# ----------------- STAGATE -----------------
Cal_Spatial_Net(adata, rad_cutoff=150)
Stats_Spatial_Net(adata)
adata = train_STAGATE(adata, alpha=0)

# ----------------- Clustering -----------------
print(f"Running mclust with k = {K}")
adata = mclust_R(adata, used_obsm='STAGATE', num_cluster=K)
adata.obs[f"mclust_{K}"] = adata.obs["mclust"].copy()
del adata.obs["mclust"]

# ----------------- Save labels -----------------
results = pd.DataFrame({
    'barcode': adata.obs_names,
    f"mclust_{K}": adata.obs[f"mclust_{K}"].astype(str).values
})
if "spatial" in adata.obsm:
    results[['x', 'y']] = adata.obsm['spatial']
else:
    results['x'] = pd.NA
    results['y'] = pd.NA

results.to_csv(OUTPUT_CSV, index=False)
print(f"✓ Saved clustering results to {OUTPUT_CSV}")

# ----------------- Compute ARI -----------------
if "ground_truth" in adata.obs.columns:
    print("\nComputing ARI scores...")
    gt = pd.Categorical(adata.obs["ground_truth"].astype(str)).codes
    pr = pd.Categorical(adata.obs[f"mclust_{K}"]).codes
    ari_val = metrics.adjusted_rand_score(gt, pr)

    print("\n" + "="*50)
    print(f"ARI (STAGATE, k={K}): {ari_val:.4f}")
    print("="*50)

    with open(OUTPUT_ARI_TXT, "w") as f:
        f.write(f"STAGATE clustering results (k={K}):\n")
        f.write("="*50 + "\n")
        f.write(f"ARI = {ari_val:.6f}\n")
    print(f"✅ Saved ARI summary: {OUTPUT_ARI_TXT}")
else:
    print("\n⚠️  No ground_truth column found — ARI computation skipped.")

print("\nDone.")
logfile.close()