import pandas as pd
import networkx as nx
from pyvis.network import Network
from collections import defaultdict
from sklearn.metrics.cluster import normalized_mutual_info_score

# Load the data
file_path = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_catElmo/results/tenClusters/tutorial--RotationEncodingBL62.txt'
data = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=['Sequence', 'ClusterID', 'Vgene', 'Info'])

# Create a graph
G = nx.Graph()

# Add nodes and edges
for _, row in data.iterrows():
    sequence = row['Sequence']
    cluster_id = row['ClusterID']
    G.add_node(sequence, title=f"Cluster: {cluster_id}, Vgene: {row['Vgene']}, Info: {row['Info']}")
    # Connect sequences in the same cluster
    for _, row2 in data[data['ClusterID'] == cluster_id].iterrows():
        if sequence != row2['Sequence']:
            G.add_edge(sequence, row2['Sequence'])

# Calculate number of clusters
clusters = list(nx.connected_components(G))
num_clusters = len(clusters)

# Calculate retention and purity fraction
pure_clusters = 0
pure_cluster_retention = 0
total_sequences = len(data)

for cluster in clusters:
    cluster_sequences = [data.loc[data['Sequence'] == seq] for seq in cluster]
    epitopes = [str(seq['Info'].values[0]).split('|')[0] for seq in cluster_sequences]
    most_common_epitope = max(set(epitopes), key=epitopes.count)
    purity = epitopes.count(most_common_epitope) / len(epitopes)
    if purity == 1.0:
        pure_clusters += 1
        pure_cluster_retention += len(cluster)

pure_cluster_fraction = pure_clusters / num_clusters
retention = pure_cluster_retention / total_sequences

# Print the analysis results
print(f"Number of Clusters: {num_clusters}")
print(f"Pure Cluster Fraction: {pure_cluster_fraction:.2%}")
print(f"Retention: {retention:.2%}")

# Create a PyVis network
net = Network(notebook=True)
net.from_nx(G)
net.show('GIANAtree_visualization.html')