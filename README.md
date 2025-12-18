# Temporal analysis subsets

Targeted treatment aims to make patients more similar to each, reducing the ability to identify subsets. We developed a Within-Subset Coherency (WSC) to assess temporal stability. This metric examines whether patients consistently cluster more closely with their baseline "neighbours"—those with similar joint involvement patterns at baseline (k=50 nearest neighbors across different clinical signatures (i.e. patietn embedding)).

## Methodology 
Our temporal analysis function works as follows"

On baseline:
-	Step 1: Run k-nearest neighbors (k = 50) to identify each patient’s baseline neighbors (according to cosine distance)

Then for each visit: 
-	Step 2: Compute the within-neighborhood distance and the between-neighborhood distance, and calculate their ratio (we call this : Within Subset Coherency).
-	Step 3: Check for every patient if we maintain a smaller ‘within’ distance with baseline neighbours compared to all other patients (ratio > 1 = good, ratio < 1 = bad)

