import pandas as pd

# Load the data
file_path = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_catElmo/TCRantigenData_unique.txt'
data = pd.read_csv(file_path, sep='\t', header=None)

# Extract the antigen column
data.columns = ['aminoAcid', 'TRBV', 'column3', 'column4', 'Antigen']
antigens = data['Antigen']

# Get the unique antigens
unique_antigens = antigens.unique()

# Output the number of unique antigens (clusters)
num_clusters = len(unique_antigens)
print(f'The number of clusters required as the ground truth is: {num_clusters}')