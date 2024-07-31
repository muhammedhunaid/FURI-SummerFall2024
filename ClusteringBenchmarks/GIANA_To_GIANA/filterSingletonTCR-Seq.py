import pandas as pd

# Define the file paths
input_file_path = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_GIANA/TCRantigenData_unique.txt'
output_file_path = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_GIANA/filteredTCRantigenData_unique.txt'

# Read the file and parse the data
with open(input_file_path, 'r') as f:
    lines = f.readlines()

# Step 2: Count the occurrences of each antigen
antigen_count = {}
for line in lines:
    antigen = line.strip().split('\t')[-1]  # Antigen is in the last column
    if antigen not in antigen_count:
        antigen_count[antigen] = 0
    antigen_count[antigen] += 1

# Step 3: Filter out Singleton TCRs (those with antigens appearing only once)
filtered_lines = []
for line in lines:
    antigen = line.strip().split('\t')[-1]
    if antigen_count[antigen] > 1:
        filtered_lines.append(line)

# Step 4: Write the output file with the filtered data
with open(output_file_path, 'w') as f:
    f.writelines(filtered_lines)

output_file_path