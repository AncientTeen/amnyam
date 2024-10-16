from __future__ import annotations

import numpy as np
from icecream import ic


def distance_matrix(buff_w: list[list[float]], metric: int) -> list[list[float]]:
    # buff_w = np.column_stack(buff)
    # ic(buff_w)
    N = len(buff_w)
    dist_matrix = None

    if metric == 1:
        # Euclidean distance
        dist_matrix = [[np.sqrt(np.sum((buff_w[j] - buff_w[i]) ** 2)) for i in range(N)] for j in range(N)]

    elif metric == 2:
        # Manhattan distance
        dist_matrix = [[np.sum(np.abs(buff_w[j] - buff_w[i])) for i in range(N)] for j in range(N)]

    elif metric == 3:
        # Chebyshev distance
        dist_matrix = [[np.max(np.abs(buff_w[j] - buff_w[i])) for i in range(N)] for j in range(N)]

    elif metric == 4:
        # Mahalanobis distance
        cov_matrix = np.cov(np.array(buff_w).T)
        inv_cov_matrix = np.linalg.inv(cov_matrix)

        dist_matrix = [
            [
                np.sqrt((buff_w[j] - buff_w[i]).T @ inv_cov_matrix @ (buff_w[j] - buff_w[i]))
                for i in range(N)
            ]
            for j in range(N)
        ]

    # Check the shape of the distance matrix
    ic(np.shape(dist_matrix))
    return dist_matrix


def lance_williams_single_linkage(dist_matrix, i, j, k):
    d_ik = dist_matrix[i][k]
    d_jk = dist_matrix[j][k]
    updated_distance = 0.5 * d_ik + 0.5 * d_jk - 0.5 * abs(d_ik - d_jk)
    return updated_distance

def single_linkage_clustering(buff, metric, num_clusters=2):
    # Step 1: Calculate the distance matrix
    buff_w = np.column_stack(buff)

    dist_matrix = distance_matrix(buff_w, metric)
    N = len(buff_w)  # Number of points

    # Initialize each point as its own cluster
    clusters = {i: [i] for i in range(N)}

    # Debug: Check the structure of clusters and distance matrix
    # ic(N, clusters, dist_matrix)

    while len(clusters) > num_clusters:
        min_dist = float('inf')
        pair_to_merge = None

        # Find the closest pair of clusters
        for i in clusters:
            for j in clusters:
                if i != j:
                    distance = dist_matrix[i][j]  # Ensure i and j are within bounds
                    if distance < min_dist:
                        min_dist = distance
                        pair_to_merge = (i, j)

        if pair_to_merge:
            i, j = pair_to_merge

            # Merge cluster j into cluster i
            clusters[i].extend(clusters[j])
            del clusters[j]

            # Update the distance matrix using Lance-Williams formula
            for k in clusters:
                if k != i:
                    dist_matrix[i][k] = lance_williams_single_linkage(dist_matrix, i, j, k)
                    dist_matrix[k][i] = dist_matrix[i][k]

            # Set distances for the merged cluster j to infinity
            for k in range(N):
                dist_matrix[j][k] = dist_matrix[k][j] = float('inf')

    # Assign points to clusters
    cluster_labels = [-1] * N
    for cluster_id, cluster_points in clusters.items():
        for point in cluster_points:
            cluster_labels[point] = cluster_id

    return cluster_labels
