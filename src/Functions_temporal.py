import numpy as np
import pandas as pd

def one_vs_rest_silhouette(D_visit, cluster_labels):
    """
    D_visit: precomputed distance matrix (n_samples x n_samples)
    cluster_labels: array/series of cluster assignments
    """
    n_samples = D_visit.shape[0]
    ovr_scores = np.zeros(n_samples)
    unique_clusters = np.unique(cluster_labels)
    
    for i in range(n_samples):
        current_cluster = cluster_labels[i]
        
        # Get indices for SAME cluster and OTHER clusters
        same_cluster_idx = np.where(cluster_labels == current_cluster)[0]
        other_clusters_idx = np.where(cluster_labels != current_cluster)[0]
        
        # --- 1. a_i (Cohesion) ---
        # Mean distance to others in same cluster (exclude self)
        mask_a = same_cluster_idx[same_cluster_idx != i]
        a_i = D_visit[i, mask_a].mean() if len(mask_a) > 0 else 0
        
        # --- 2. b_i (One-vs-Rest Separation) ---
        # Mean distance to EVERYONE not in the current cluster
        if len(other_clusters_idx) > 0:
            b_i = D_visit[i, other_clusters_idx].mean()
        else:
            b_i = 0
            
        # --- 3. s_i (One-vs-Rest Score) ---
        ovr_scores[i] = (b_i - a_i) / max(a_i, b_i) if max(a_i, b_i) > 0 else 0
        
    return ovr_scores

def manual_silhouette_samples(D_visit, cluster_labels):
    """
    D_visit: precomputed distance matrix (n_samples x n_samples)
    cluster_labels: array/series of cluster assignments
    """
    n_samples = D_visit.shape[0]
    individual_scores = np.zeros(n_samples)
    unique_clusters = np.unique(cluster_labels)
    
    # Create a mapping of cluster names to their indices for speed
    cluster_indices = {c: np.where(cluster_labels == c)[0] for c in unique_clusters}

    for i in range(n_samples):
        current_cluster = cluster_labels[i]
        members_same_cluster = cluster_indices[current_cluster]
        
        # --- 1. Calculate a_i ---
        # Mask out the point itself (distance to self is 0)
        mask_a = members_same_cluster[members_same_cluster != i]
        
        if len(mask_a) > 0:
            a_i = D_visit[i, mask_a].mean()
        else:
            # Silhouette is technically 0 if a cluster has only 1 member
            individual_scores[i] = 0
            continue

        # --- 2. Calculate b_i ---
        other_cluster_averages = []
        for other_cluster in unique_clusters:
            if other_cluster == current_cluster:
                continue
                
            members_other_cluster = cluster_indices[other_cluster]
            if len(members_other_cluster) > 0:
                dist_to_other = D_visit[i, members_other_cluster].mean()
                other_cluster_averages.append(dist_to_other)
        
        # Separation is the distance to the NEAREST neighboring cluster
        b_i = min(other_cluster_averages) if other_cluster_averages else 0

        # --- 3. Calculate s_i ---
        # Formula: (b - a) / max(a, b)
        individual_scores[i] = (b_i - a_i) / max(a_i, b_i) if max(a_i, b_i) > 0 else 0
        
    return individual_scores