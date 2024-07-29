import pandas as pd
from collections import defaultdict

# Load the clustering result file
file_path = 'tmp/TestReal-ADIRP0000023_TCRB_query_hc10s10.txt'
data = pd.read_csv(file_path, sep='\t', header=None)
data.columns = ['TCR', 'ClusterID', 'GeneSegment', 'SimilarityScore', 'Probability', 'FileSource', 'DatasetType']

def calculate_purity_retention(data):
    # Dictionary to hold clusters
    clusters = defaultdict(list)

    # Group TCRs by cluster
    for _, row in data.iterrows():
        clusters[row['ClusterID']].append((row['TCR'], row['FileSource']))

    total_tcrs = len(data)
    pure_clusters = 0
    pure_tcrs_count = 0

    for cluster_id, tcrs in clusters.items():
        source_counts = defaultdict(int)

        for _, source in tcrs:
            source_counts[source] += 1

        # Determine the most common source in the cluster
        most_common_source_count = max(source_counts.values())
        cluster_purity = most_common_source_count / len(tcrs)

        # Check if the cluster is pure
        if cluster_purity == 1:
            pure_clusters += 1
            pure_tcrs_count += len(tcrs)

    # Calculate metrics
    purity = pure_clusters / len(clusters)
    retention = pure_tcrs_count / total_tcrs

    return purity, retention

# Calculate the purity and retention
purity, retention = calculate_purity_retention(data)

print(f"Purity: {purity:.2f}")
print(f"Retention: {retention:.2f}")