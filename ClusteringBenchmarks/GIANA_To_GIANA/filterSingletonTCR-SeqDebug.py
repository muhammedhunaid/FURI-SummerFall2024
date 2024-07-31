import pandas as pd

# Define the file paths
input_file_path = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_GIANA/TCRantigenData_unique.txt'
#output_file_path_with_summary = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_GIANA/filteredTCRantigenData_unique-Debug.txt'

# Reuse the previous steps to read the input file and parse the data

# Read the file and parse the data
with open(input_file_path, 'r') as f:
    lines = f.readlines()

# Count the occurrences of each antigen
antigen_count = {}
for line in lines:
    antigen = line.strip().split('\t')[-1]  # Antigen is in the last column
    if antigen not in antigen_count:
        antigen_count[antigen] = 0
    antigen_count[antigen] += 1

# Filter out Singleton TCRs (those with antigens appearing only once)
filtered_lines = []
singleton_antigens = []
for line in lines:
    antigen = line.strip().split('\t')[-1]
    if antigen_count[antigen] > 1:
        filtered_lines.append(line)
    elif antigen_count[antigen] == 1:
        singleton_antigens.append(antigen)

# Calculate the number of unique antigens appearing at most once
unique_antigen_count = len(singleton_antigens)

# Write the output file with the filtered data and summary
output_file_path_with_singletons = '/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_GIANA/filteredTCRantigenData_unique-Debug.txt'
with open(output_file_path_with_singletons, 'w') as f:
    f.writelines(filtered_lines)
    # Append the summary and list of Singleton antigens
    f.write(f"\n\nTotal antigens appearing at most once (Singleton antigens): {unique_antigen_count}\n")
    f.write("No Singleton TCRs (antigens appearing only once) are present in this file.\n\n")
    f.write("List of Singleton antigens:\n")
    for antigen in singleton_antigens:
        f.write(f"{antigen}\n")

output_file_path_with_singletons