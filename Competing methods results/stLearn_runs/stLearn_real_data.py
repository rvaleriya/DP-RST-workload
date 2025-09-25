import stlearn as st
import scanpy as sc
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
st.settings.set_figure_params(dpi=180)

logfile = open("stLearn_RealData_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

# Define datasets and cluster numbers to try per dataset
dataset_config = {
    "Brain": [5, 7, 9],
    "Breast": [3, 5, 7],
    "Gut": [3, 5, 7, 9, 12],
    "Gut_reduced": [3, 5, 7, 9, 12],
    "Prostate": [3, 5, 7]  
}

# specify PATH to data
BASE_PATH = Path("../ST_Datasets")

# spot tile is the intermediate result of image pre-processing
TILE_PATH = Path("/tmp/tiles")
TILE_PATH.mkdir(parents=True, exist_ok=True)

# output path
OUT_PATH = Path("res_realdata_plots")
OUT_PATH.mkdir(parents=True, exist_ok=True)

# Create directory for CSV results
CSV_PATH = Path("res_realdata_csv")
CSV_PATH.mkdir(parents=True, exist_ok=True)

def run_analysis(data, n_pcs, n_clusters, use_image=True):
    # First run PCA
    st.em.run_pca(data, n_comps=n_pcs)
    
    data_processed = data.copy()
    
    if use_image:
        # Apply stSME to normalise log transformed data
        st.spatial.SME.SME_normalize(data_processed, use_data="raw")
        data_processed.X = data_processed.obsm['raw_SME_normalized']
        st.pp.scale(data_processed)
    else:
        # Use standard scaling without spatial information
        st.pp.scale(data_processed)
    
    st.em.run_pca(data_processed, n_comps=n_pcs)
    
    # K-means clustering
    st.tl.clustering.kmeans(data_processed, n_clusters=n_clusters, use_data="X_pca", key_added="kmeans")
    return data_processed

def save_clustering_results(adata, dataset_name, n_pcs, n_clusters, use_image=True):
    # Create DataFrame with spatial coordinates and cluster assignments
    results_df = pd.DataFrame({
        'spot': adata.obs_names,
        'x': adata.obsm['spatial'][:, 0],
        'y': adata.obsm['spatial'][:, 1],
        f'cluster_{n_clusters}': adata.obs['kmeans']  # Column name includes number of clusters
    })
    
    # Save to CSV
    image_suffix = "with_image" if use_image else "no_image"
    results_df.to_csv(CSV_PATH / f"{dataset_name}_clustering_{n_pcs}PCs_{n_clusters}clusters_{image_suffix}.csv", index=False)
    return results_df

def process_dataset(dataset_name, n_clusters_list, use_image=True):
    print(f"\nProcessing dataset: {dataset_name}")
    
    # Load data
    data = st.Read10X(BASE_PATH / dataset_name)
    print(data)
    
    # If Gut_reduced: subset spots based on CSV
    if dataset_name == "Gut_reduced":
        print("Subsetting Gut_reduced to selected spots...")
        subset_csv = BASE_PATH / dataset_name / "gut_df_wt_muscle_rownames.csv"
        spot_ids = pd.read_csv(subset_csv, header=None).squeeze().astype(str)
        data = data[data.obs_names.isin(spot_ids)].copy()
        print(f"✓ Subsetted to {data.n_obs} spots")
        print(data)
    
    # Pre-processing for gene count table
    st.pp.normalize_total(data)
    st.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=2000)
    data = data[:, data.var.highly_variable]
    print(data)
    
    # Pre-processing for spot image (optional)
    if use_image:
        st.pp.tiling(data, TILE_PATH)
        st.pp.extract_feature(data)
    
    # Initialize DataFrames for each PC
    results_3pcs = pd.DataFrame({
        'spot': data.obs_names,
        'x': data.obsm['spatial'][:, 0],
        'y': data.obsm['spatial'][:, 1]
    })
    results_10pcs = results_3pcs.copy()
    
    # Run analysis for different numbers of clusters
    for n_clusters in n_clusters_list:
        print(f"Running analysis with {n_clusters} clusters")
        data_3pcs = run_analysis(data, n_pcs=3, n_clusters=n_clusters, use_image=use_image)
        data_10pcs = run_analysis(data, n_pcs=10, n_clusters=n_clusters, use_image=use_image)
        
        # Save individual clustering results
        results_3pcs_clusters = save_clustering_results(data_3pcs, dataset_name, 3, n_clusters, use_image)
        results_10pcs_clusters = save_clustering_results(data_10pcs, dataset_name, 10, n_clusters, use_image)
        
        # Add clustering results to respective DataFrames
        results_3pcs[f'cluster_{n_clusters}'] = results_3pcs_clusters[f'cluster_{n_clusters}']
        results_10pcs[f'cluster_{n_clusters}'] = results_10pcs_clusters[f'cluster_{n_clusters}']
        
        # Create comparison plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        st.pl.cluster_plot(data_3pcs, use_label="kmeans", ax=ax1, size=25)
        ax1.set_title(f"3 PCs, {n_clusters} clusters")
        st.pl.cluster_plot(data_10pcs, use_label="kmeans", ax=ax2, size=25)
        ax2.set_title(f"10 PCs, {n_clusters} clusters")
        plt.suptitle(f"{dataset_name} Dataset")
        plt.tight_layout()
        
        # Add image info to filename
        image_suffix = "with_image" if use_image else "no_image"
        plt.savefig(str(OUT_PATH / f"{dataset_name}_clustering_{n_clusters}clusters_comparison_{image_suffix}.png"), 
                   dpi=180, bbox_inches='tight')
        plt.close()
    
    # Save all results for each PC in a single file
    image_suffix = "with_image" if use_image else "no_image"
    results_3pcs.to_csv(CSV_PATH / f"{dataset_name}_clustering_3PCs_all_clusters_{image_suffix}.csv", index=False)
    results_10pcs.to_csv(CSV_PATH / f"{dataset_name}_clustering_10PCs_all_clusters_{image_suffix}.csv", index=False)

# Process each dataset
for dataset_name, n_clusters_list in dataset_config.items():
    # Run with image information
    process_dataset(dataset_name, n_clusters_list, use_image=True)
    # Run without image information
    process_dataset(dataset_name, n_clusters_list, use_image=False)