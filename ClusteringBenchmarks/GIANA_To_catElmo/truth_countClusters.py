import pandas as pd

# Load the data from the file
file_path = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_catElmo/tutorial.txt'
data = pd.read_csv(file_path, sep='\t')

# Extract the unique antigens (TRBV values)
unique_antigens = data['vMaxResolved'].unique()

# Output the number of unique antigens
num_clusters = len(unique_antigens)
print(f'The number of clusters required as the ground truth is: {num_clusters}')