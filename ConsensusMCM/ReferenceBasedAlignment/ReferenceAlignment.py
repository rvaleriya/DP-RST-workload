import numpy as np
import pandas as pd
import os
import glob
import re
import argparse
from scipy.linalg import sqrtm

def gaussian_wasserstein_distance(mean1, cov1, mean2, cov2, sqrt_cov2=None):
    """Computes the 2-Wasserstein distance between two multivariate Gaussians."""
    mean_diff = mean1 - mean2
    mean_term = np.dot(mean_diff, mean_diff)

    if sqrt_cov2 is None:
        sqrt_cov2 = sqrtm(cov2)

    if np.iscomplexobj(sqrt_cov2):
        sqrt_cov2 = np.real_if_close(sqrt_cov2)

    inner = sqrtm(sqrt_cov2 @ cov1 @ sqrt_cov2)

    if np.iscomplexobj(inner):
        inner = np.real_if_close(inner)

    cov_term = np.trace(cov1 + cov2 - 2.0 * inner)
    cov_term = np.clip(cov_term, a_min=0.0, a_max=None)

    distance_sq = mean_term + cov_term
    return np.sqrt(max(distance_sq, 0.0))

def load_real_data(base_dir):
    """Loads data from DP-RST post-processing output directories."""
    print(f"--- Loading Real Data from: {base_dir} ---")
    shard_results = []
    
    means_dir = os.path.join(base_dir, 'means')
    cov_dir = os.path.join(base_dir, 'covariance')
    prop_dir = os.path.join(base_dir, 'proportions')

    if not os.path.isdir(prop_dir):
        raise FileNotFoundError(f"Proportions directory not found: {prop_dir}")

    prop_files = sorted(glob.glob(os.path.join(prop_dir, "Shard_*_proportions.csv")))
    if not prop_files:
        raise FileNotFoundError(f"No proportion files found in {prop_dir}")

    for prop_path in prop_files:
        shard_id_match = re.search(r"Shard_(\d+)_proportions\.csv", os.path.basename(prop_path))
        if not shard_id_match:
            continue
        shard_id_str = shard_id_match.group(1)

        means_path = os.path.join(means_dir, f"Shard_{shard_id_str}_means.csv")
        cov_path = os.path.join(cov_dir, f"Shard_{shard_id_str}_covariance.csv")

        if not os.path.exists(means_path) or not os.path.exists(cov_path):
            print(f"Warning: Missing means or covariance file for shard {shard_id_str}. Skipping.")
            continue
            
        try:
            prop_df = pd.read_csv(prop_path)
            means_df = pd.read_csv(means_path)
            cov_matrix = pd.read_csv(cov_path).values

            merged_df = pd.merge(prop_df, means_df, left_on='cluster', right_on='team')

            if merged_df.empty:
                print(f"Warning: No matching clusters for shard {shard_id_str}. Skipping.")
                continue

            n_obs = merged_df['count'].values
            mean_vectors = merged_df.drop(columns=['cluster', 'count', 'proportion', 'total_n', 'team']).values
            
            n_local_clusters = len(n_obs)
            covariances = np.tile(cov_matrix[np.newaxis, :, :], (n_local_clusters, 1, 1))

            shard_results.append({
                'shard_id': int(shard_id_str),
                'means': mean_vectors,
                'covariances': covariances,
                'n_obs': n_obs
            })
        except Exception as e:
            print(f"Error processing shard {shard_id_str}: {e}")
            
    if not shard_results:
        raise ValueError("No shard data could be loaded. Please check file paths and formats.")

    n_features = shard_results[0]['means'].shape[1]
    print(f"Successfully loaded data for {len(shard_results)} shards.")
    print(f"Inferred number of features (PCs): {n_features}")
    
    return shard_results

def align_to_reference(shard_results, ref_idx):
    """
    Aligns all local clusters to a reference shard using the nearest
    centroid based on Wasserstein distance.
    """
    reference_shard = shard_results[ref_idx]
    global_means = reference_shard['means']
    global_covs = reference_shard['covariances']
    n_global = len(global_means)

    global_cov_sqrts = [
        np.real_if_close(sqrtm(global_covs[k])) for k in range(n_global)
    ]
    
    print(f"--- Aligning to Reference Shard {reference_shard['shard_id']} ({n_global} clusters) ---")

    all_hard_mappings = []
    all_soft_probs = []
    all_local_means = []
    
    for shard in shard_results:
        local_centroids = shard['means']
        local_covariances = shard['covariances']
        
        cost_matrix = np.zeros((local_centroids.shape[0], n_global))

        for j, (local_mean, local_cov) in enumerate(zip(local_centroids, local_covariances)):
            for k in range(n_global):
                cost_matrix[j, k] = gaussian_wasserstein_distance(
                    local_mean, local_cov,
                    global_means[k], global_covs[k],
                    sqrt_cov2=global_cov_sqrts[k]
                )
        
        # Convert distances to probabilities via softmax on negative cost
        temperature = 1.0
        stabilized_cost = cost_matrix - cost_matrix.min(axis=1, keepdims=True)
        soft_probs = np.exp(-stabilized_cost / temperature)
        soft_probs /= soft_probs.sum(axis=1, keepdims=True)
        
        hard_mapping = np.argmax(soft_probs, axis=1)
        
        all_hard_mappings.append(hard_mapping)
        all_soft_probs.append(soft_probs)
        all_local_means.append(local_centroids)
        
    hard_assignments = np.concatenate(all_hard_mappings)
    soft_probabilities = np.vstack(all_soft_probs)
    all_local_means = np.vstack(all_local_means)

    avg_confidence = soft_probabilities.max(axis=1).mean()
    print(f"Average assignment confidence: {avg_confidence:.3f}")
    
    # Calculate global weights and update means based on assignments
    all_n_obs = np.concatenate([s['n_obs'] for s in shard_results])
    global_weights = np.zeros(n_global)
    for k in range(n_global):
        mask = hard_assignments == k
        if mask.any():
            global_weights[k] = all_n_obs[mask].sum()
            global_means[k] = np.average(
                all_local_means[mask],
                weights=all_n_obs[mask],
                axis=0
            )
    global_weights /= global_weights.sum()

    global_params = {
        'means': global_means,
        'covariances': global_covs, # Note: Covariances are from the reference shard
        'weights': global_weights
    }
    
    return hard_assignments, soft_probabilities, global_params


def main(args):
    """Main execution function."""
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"All outputs will be saved to: {args.output_dir}")
    
    shard_data = load_real_data(args.data_dir)
    
    # Find the index of the user-specified reference shard
    shard_ids = [s['shard_id'] for s in shard_data]
    try:
        ref_idx = shard_ids.index(args.ref_shard)
    except ValueError:
        raise ValueError(f"Reference shard ID {args.ref_shard} not found. Available IDs: {shard_ids}")
    
    hard_assignments, soft_probabilities, global_params = align_to_reference(
        shard_data, ref_idx=ref_idx
    )
    
    # Save global parameters
    print("\n--- Saving Global Cluster Parameters ---")
    
    np.save(os.path.join(args.output_dir, 'global_weights.npy'), global_params['weights'])
    np.save(os.path.join(args.output_dir, 'global_means.npy'), global_params['means'])
    np.save(os.path.join(args.output_dir, 'global_covariances.npy'), global_params['covariances'])
    
    print(f"Global parameters saved to: {args.output_dir}")
    print("Weights (rounded):\n", np.round(global_params['weights'], 3))
    
    # Save mappings
    print("\n--- Saving Local-to-Global Mapping ---")
    
    hard_assign_path = os.path.join(args.output_dir, 'local_to_global_hard_assignments.csv')
    soft_assign_path = os.path.join(args.output_dir, 'local_to_global_soft_probabilities.csv')
    
    pd.DataFrame(
        hard_assignments, columns=['global_cluster_assignment']
    ).to_csv(hard_assign_path, index_label='local_cluster_index')
    
    pd.DataFrame(soft_probabilities).to_csv(
        soft_assign_path, index_label='local_cluster_index'
    )

    print(f"Mapping files saved to: {args.output_dir}")
    print("\n--- Alignment Complete ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reference-based alignment for metaclustering using Wasserstein distance."
    )
    parser.add_argument(
        '--data-dir', 
        type=str,
        required=True,
        help='Directory containing the post-processed DP-RST results.'
    )
    parser.add_argument(
        '--output-dir', 
        type=str, 
        required=True,
        help='Directory to save the alignment results.'
    )
    parser.add_argument(
        '--ref-shard',
        type=int,
        required=True,
        help='The ID of the shard to use as the reference for alignment.'
    )
    
    args = parser.parse_args()
    main(args)
