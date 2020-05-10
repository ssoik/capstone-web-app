import numpy as np
import pandas as pd
from sklearn.cluster import AffinityPropagation as ap

# Compute a similarity matrix

def similarity_matrix(id_series, filepath):
    # Rename for merging if necessary
    if id_series.name != 'drugbank-id':
        id_series = pd.Series(id_series, name='drugbank-id')

    # Compute matrix
    matrix = []
    for i, drug in enumerate(id_series):
        inds = pd.read_csv(f'{filepath}sim-data/{drug}.txt')
        matrix.append(pd.merge(id_series, inds, on='drugbank-id')['dice-similarity'])
        
    return np.array(matrix)

# Find average similarity by cluster using Affinity Propagation
# Set similarity threshold from upper half of cluster averages

def similarity_threshold(sim_matrix, id_series, filepath, cutoff=0.5):
    # If zero or one drug found, return cutoff
    if len(sim_matrix) <= 1:
        return cutoff
    else:
        # Find clusters
        clusters = ap().fit(sim_matrix)
        cluster_centers_indices = clusters.cluster_centers_indices_
        labels = clusters.labels_
        n_clusters = len(cluster_centers_indices)

        # Compute mean similarity for each cluster
        k_means = []
        for k in range(n_clusters):
            class_members = labels == k
            member_ids = id_series[class_members]
            member_matrix = similarity_matrix(member_ids, filepath)

            if len(member_matrix) == 1:
                k_means.append(np.mean(member_matrix))
            else:
                uptri = []
                for i, row in enumerate(member_matrix):
                    uptri += list(row[i + 1:])

                k_mean = np.mean(uptri)
                k_means.append(k_mean)

        # Take mean of upper half of clusters
        thresh = np.mean(sorted(k_means)[len(k_means) // 2:])
        if thresh >= cutoff and thresh < 1:
            return thresh
        else:
            return cutoff
