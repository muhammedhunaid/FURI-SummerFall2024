import pandas as pd
import networkx as nx
from pyvis.network import Network

# Load the data
file_path = '/Users/MuhammedH/github/GIANA/results/exactMode/TRBVgeneAvailable/tutorial--RotationEncodingBL62.txt'
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

# Create a PyVis network
net = Network(notebook=True)
net.from_nx(G)
net.show('GIANAtree_visualization.html')