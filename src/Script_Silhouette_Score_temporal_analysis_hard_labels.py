from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import pairwise_distances
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import silhouette_samples

k = 50  # define neighbourhood size
CLUSTER = 'Original_cluster' # = list with all clinical signatures # (i.e. lymphoid/ myeloid levels) - which indicates baseline profile of patient.
PATIENT_ID = 'PATNR' # variable that defines patients
VISIT_ID = 'visnr'                
    
# --- Hard Label clustering: Define neighbors by Baseline Cluster instead of KNN ---
# Create a dictionary where keys are cluster IDs and values are lists of PATNRs
cluster_map = df_raw[df_raw[VISIT_ID] == 1].set_index('PATNR')[CLUSTER].to_dict()

for time in sorted(df_raw[VISIT_ID].unique()):
    if time > 4: continue

    # Get data for this visit and map their original clusters
    df_visit = df_raw[df_raw[VISIT_ID] == time].copy().set_index(PATIENT_ID)
    df_visit[CLUSTER] = df_visit.index.map(cluster_map)
    
    # Drop patients who weren't there at baseline or have no cluster
    df_visit = df_visit.dropna(subset=[CLUSTER])
    
    if len(df_visit) < 2: continue

    # Calculate full distance matrix for the visit
    X_visit = df_visit[lf_col].values
    D_visit = pairwise_distances(X_visit, metric='cosine')
    
    # Get indices for each cluster present at this visit
    cluster_indices = df_visit.groupby(CLUSTER).indices
    
    # Run the silhouette-like function (cluster vs rest)
    individual_ss = ft.one_vs_rest_silhouette(D_visit, df_visit[CLUSTER].values)
    # or true silhouette (cluster vs nearest cluster)
    #individual_ss = ft.manual_silhouette_samples(D_visit, df_visit[CLUSTER].values)

    # 2. Add these scores back to your visit dataframe so we can group them
    df_visit['patient_ss'] = individual_ss

    # 3. Now, when you iterate through your clusters:
    for cluster_id, member_indices in cluster_indices.items():
        if len(member_indices) < 2: continue
            
        # Get the individual scores for just this cluster
        cluster_scores = df_visit[df_visit[CLUSTER] == cluster_id]['patient_ss']
        
        for score in cluster_scores:
            results.append({
                "Time Point": time,
                "Cluster": cluster_id,
                "Individual_SS": score  # Raw score per patient
            })

            
# -------------------
# Plot temporal stability -> assuming four subgroups
# -------------------
plt.figure()
pal_colorblind = ['#4F6CCF', '#2db9cc', '#FBC93D', '#FA4D4D'];
# Now Seaborn works perfectly:
sns.lineplot(data=pd.DataFrame(results), x='Time Point', y='Individual_SS', 
             hue='Cluster', errorbar='sd', marker='o', linewidth=4, palette=pal_colorblind)

# Optional: Customize the plot
# Example time points = 1,2,3,4 
therapy = 'y'
plt.xticks([1, 2, 3, 4], ["0", "3", "6", "12"], fontsize=14) 
plt.yticks(fontsize=14)
plt.title("Silhouette Score per cluster", fontsize=14)
plt.xlabel("Time (months)", fontsize=14) 
plt.ylabel('Silhouette Score', fontsize=14)  #  (Ratio BCSS vs WCSS)
plt.legend(title='Hard label cluster',fontsize=14) # , []
plt.axhline(0, color='gray', linestyle='--', linewidth=3)
plt.ylim([-.5,1])
plt.grid(True)

pal_colorblind = ['#4F6CCF', '#2db9cc', '#FBC93D', '#FA4D4D'];
df_stability = pd.DataFrame(results)
