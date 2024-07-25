import pandas as pd

# Load the data
file_path = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_catElmo/trimmedShantanuAntigenData.txt'
data = pd.read_csv(file_path, sep='\t', header=None)

# Assign column names based on the previous analysis
data.columns = ['aminoAcid', 'TRBV', 'column3', 'column4', 'Antigen']

# Filter sequences with length 15 characters
filtered_data = data[data['aminoAcid'].str.len() == 15]

# Count unique antigens for these sequences
unique_antigens = filtered_data['Antigen'].unique()
num_unique_antigens = len(unique_antigens)

print(f'The number of unique antigens for TCR sequences with length 15 characters is: {num_unique_antigens}')