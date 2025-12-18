from sklearn.neighbors import NearestNeighbors
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist, squareform 
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from sklearn.metrics.pairwise import cosine_distances, pairwise_distances
from scipy import stats

k = 50  # define neighbourhood size
lf_col # = list with all clinical signatures # (i.e. lymphoid/ myeloid levels)

# Import dataframe with Latent factors (clinical signatures) + clusters from NORDSTAR
df_raw  = df_meta_rep.copy()


# First define baseline neighbours (take subset of visit 1)
df_base = df_raw[df_raw["visnr"] == 1].copy().set_index("PATNR")
X_base = df_base[lf_col].values

# fit kNN
n_neighbors = min(k+1, len(df_base))   # +1 includes self
knn = NearestNeighbors(n_neighbors=n_neighbors, metric='cosine')
knn.fit(X_base)

dist_base, neigh_base = knn.kneighbors(X_base)

# remove "self" (first entry)
dist_base = dist_base[:, 1:]
neigh_base = neigh_base[:, 1:]

# mapping
row_to_pid = dict(enumerate(df_base.index))
pid_to_row = {pid: i for i, pid in enumerate(df_base.index)}

# store k-neighbour list per patient
baseline_neighbours = {}

for pid in df_base.index:
    i0 = pid_to_row[pid]
    nbr_ids = [row_to_pid[j] for j in neigh_base[i0]]
    baseline_neighbours[pid] = nbr_ids

results = []

for time in df_raw["visnr"].unique():
    print(time)
    if time in [0, 5.1]:   # skip baseline & special measurement
        continue

    df_visit = df_raw[df_raw["visnr"] == time].copy().set_index("PATNR")

    if len(df_visit) < 2:
        continue

    X_visit = df_visit[lf_col].values
    ids_visit = df_visit.index.values

    # distance matrix for this visit
    D_visit = pairwise_distances(X_visit, metric='cosine')

    # mapping for this visit
    row_visit = {pid: i for i, pid in enumerate(ids_visit)}

    for pid in ids_visit:

        # skip if not present at baseline
        if pid not in baseline_neighbours:
            continue
        
        i = row_visit[pid]
        cluster_i = df_visit.loc[pid, "Original_cluster"]

        # baseline neighbour list
        base_nbrs = baseline_neighbours[pid]

        # neighbours present at this visit
        nbrs_present = [p for p in base_nbrs if p in row_visit]

        # all non-neighbours at baseline
        non_nbrs_present = [p for p in ids_visit if p not in nbrs_present and p != pid]

        # -------------------
        # Within = distance to baseline neighbours
        # -------------------
        if len(nbrs_present) > 0:
            idx_within = [row_visit[p] for p in nbrs_present]
            dist_within = D_visit[i, idx_within].mean()
        else:
            dist_within = np.nan

        # -------------------
        # Between = distance to baseline NON-neighbours
        # -------------------
        if len(non_nbrs_present) > 0:
            idx_between = [row_visit[p] for p in non_nbrs_present]
            dist_between = D_visit[i, idx_between].mean()
        else:
            dist_between = np.nan

        # -------------------
        # Ratio
        # -------------------
        ratio = (
            dist_between / dist_within 
            if dist_within > 0 and dist_between > 0
            else np.nan
        )

        results.append({
            "patient_id": pid,
            "Cluster": cluster_i,
            "Time Point": time,
            "Within_distance": dist_within,
            "Between_distance": dist_between,
            "Ratio_between_within": ratio,
            "N_within_present": len(nbrs_present),
            "N_between_present": len(non_nbrs_present)
        })

pal_colorblind = ['#4F6CCF', '#2db9cc', '#FBC93D', '#FA4D4D'];
df_stability = pd.DataFrame(results)

sns.lineplot(x='Time Point', y='Ratio_between_within', hue='Cluster', data=df_stability, marker='o', linewidth=2.5, palette=pal_colorblind)

# Optional: Customize the plot
therapy = 'y'
plt.xticks([1, 2, 3, 4], ["0", "3", "6", "12"], fontsize=14) 
plt.yticks([0, 1, 2, 3, 4], fontsize=14)
plt.title("Group Coherence index per JIP in the patient embedding (k=50)", fontsize=14)
plt.xlabel("Time (months)", fontsize=14) 
plt.ylabel('Group coherence index', fontsize=14)  #  (Ratio BCSS vs WCSS)
plt.legend(fontsize=14)
plt.axhline(1, color='gray', linestyle='--', linewidth=3)
#plt.ylim([0,4])
plt.grid(True)
plt.legend(title='Group', fontsize=9)