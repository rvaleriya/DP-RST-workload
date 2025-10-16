# ======== stLearn on STARmap with SME normalization ========
import numpy as np
import pandas as pd
from pathlib import Path
import scanpy as sc
import stlearn as st
from sklearn import metrics
import warnings
warnings.filterwarnings("ignore")

# ---------- Paths ----------
H5AD = "/scratch/user/varogovchenko/BASTION_HPRC/STARmap/STARmap_20180505_BY3_1k_20251008011714.h5ad"
OUT  = Path("res_starmap"); OUT.mkdir(exist_ok=True)

K = 7
PCS_LIST = [3, 10]
N_COMPS = max(PCS_LIST)
np.random.seed(42)

# ---------- Load data ----------
print("Loading data...")
adata = sc.read_h5ad(H5AD)

# Ensure spatial coordinates are in obsm
if "spatial" not in adata.obsm:
    adata.obsm["spatial"] = adata.obs[["x","y"]].to_numpy()

# Store original counts
adata.layers["counts"] = adata.X.copy()

# Create dummy Visium-like fields (required by stLearn)
xy = adata.obsm["spatial"]
adata.obs["imagerow"] = np.round(xy[:,1]).astype(int)
adata.obs["imagecol"] = np.round(xy[:,0]).astype(int)
adata.obs["array_row"] = adata.obs["imagerow"]
adata.obs["array_col"] = adata.obs["imagecol"]

# Create dummy morphology features (required by stLearn, even if not used)
if "X_morphology" not in adata.obsm:
    adata.obsm["X_morphology"] = np.zeros((adata.n_obs, 8), dtype=float)

print(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")

# ---------- Pre-processing (following stLearn tutorial) ----------
print("\nPre-processing gene count table...")

# Filter, normalize, log transform (as in tutorial)
st.pp.filter_genes(adata, min_cells=1)
st.pp.normalize_total(adata)
st.pp.log1p(adata)

print("Pre-processing complete")

# ---------- Run initial PCA (needed by SME) ----------
print(f"\nRunning initial PCA with {N_COMPS} components...")
st.em.run_pca(adata, n_comps=N_COMPS)
print("Initial PCA complete")

# ---------- Apply SME normalization ----------
print("\nApplying stSME normalization...")

# Create a copy for SME processing (as in tutorial)
adata_sme = adata.copy()

# Apply SME normalize on the log-transformed data
st.spatial.SME.SME_normalize(adata_sme, use_data="raw", platform="Visium", weights="weights_matrix_pd_gd")

# Guard against NaNs/Infs from SME
adata_sme.obsm["raw_SME_normalized"] = np.nan_to_num(
    adata_sme.obsm["raw_SME_normalized"], 
    nan=0.0, 
    posinf=0.0, 
    neginf=0.0
)

# Set SME-normalized data as main matrix (as in tutorial)
adata_sme.X = adata_sme.obsm['raw_SME_normalized']

print("stSME normalization complete")

# ---------- Scale and run PCA on SME-normalized data ----------
print("\nScaling SME-normalized data...")
st.pp.scale(adata_sme)

# Ensure no NaNs/Infs after scaling
adata_sme.X = np.nan_to_num(adata_sme.X, nan=0.0, posinf=0.0, neginf=0.0)

print(f"Running final PCA with {N_COMPS} components on SME-normalized data...")
st.em.run_pca(adata_sme, n_comps=N_COMPS)
print("Final PCA complete")

# ---------- Run clustering for different PC counts ----------
print("\nRunning k-means clustering...")

for n_pcs in PCS_LIST:
    print(f"  Clustering with {n_pcs} PCs, k={K}...")
    
    # Extract subset of PCs for this run
    adata_sme.obsm["X_pca_subset"] = adata_sme.obsm["X_pca"][:, :n_pcs].copy()
    
    # Run k-means clustering (as in tutorial)
    st.tl.clustering.kmeans(
        adata_sme, 
        n_clusters=K, 
        use_data="X_pca_subset", 
        key_added=f"label_{n_pcs}PCs"
    )
    
    print(f"    Clustering complete for {n_pcs} PCs")

# ---------- Save combined results ----------
print("\nSaving results...")
xy = adata_sme.obsm["spatial"]
out = pd.DataFrame({
    "ID": adata_sme.obs_names,
    "x": xy[:,0],
    "y": xy[:,1],
    f"label_3PCs_k{K}": adata_sme.obs["label_3PCs"],
    f"label_10PCs_k{K}": adata_sme.obs["label_10PCs"]
})

out_path = OUT / f"stlearn_STARmap_3PCs_10PCs.csv"
out.to_csv(out_path, index=False)
print(f"✅ Saved combined labels: {out_path}")

# ---------- Compute and save ARI ----------
if "ground_truth" in adata_sme.obs.columns:
    print("\nComputing ARI scores...")
    
    gt = pd.Categorical(adata_sme.obs["ground_truth"].astype(str)).codes
    pr3 = pd.Categorical(adata_sme.obs["label_3PCs"]).codes
    pr10 = pd.Categorical(adata_sme.obs["label_10PCs"]).codes
    
    ari3 = metrics.adjusted_rand_score(gt, pr3)
    ari10 = metrics.adjusted_rand_score(gt, pr10)
    
    print(f"\n{'='*50}")
    print(f"RESULTS:")
    print(f"{'='*50}")
    print(f"  ARI (3 PCs, k={K}):  {ari3:.4f}")
    print(f"  ARI (10 PCs, k={K}): {ari10:.4f}")
    print(f"{'='*50}")
    
    # Save ARI summary
    ari_path = OUT / "stlearn_STARmap_ARI.txt"
    with open(ari_path, "w") as f:
        f.write(f"stSME-based clustering results (k={K}):\n")
        f.write(f"="*50 + "\n")
        f.write(f"3PCs:  ARI = {ari3:.6f}\n")
        f.write(f"10PCs: ARI = {ari10:.6f}\n")
    
    print(f"✅ Saved ARI summary: {ari_path}")
else:
    print("\n⚠️  No ground_truth column found — ARI computation skipped.")

# ---------- Save the processed AnnData object ----------
# adata_path = OUT / "stlearn_STARmap_SME_processed.h5ad"
# adata_sme.write(adata_path)
# print(f"✅ Saved processed AnnData: {adata_path}")

print("\n stLearn clustering pipeline complete!")