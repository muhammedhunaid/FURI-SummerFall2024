import csv
from collections import defaultdict, Counter

def calculate_purity_and_retention(clustered_file_path, reference_file_path):
    clusters = defaultdict(list)
    total_tcrs = 0
    total_reference_tcrs = 0

    # Read the reference file to count the total number of TCR sequences
    with open(reference_file_path, 'r') as ref_file:
        for _ in ref_file:
            total_reference_tcrs += 1

    # Read the clustered file and parse clusters
    with open(clustered_file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('##'):
                continue  # Skip metadata lines
            cdr3_seq = row[0]
            cluster_id = row[1]
            antigen = row[-1].split(':')[-1]
            clusters[cluster_id].append((cdr3_seq, antigen))
            total_tcrs += 1

    sum_most_common = 0
    cluster_purity_dict = {}

    # Calculate purity and retention and output the results
    for cluster_id, tcr_list in clusters.items(): #what does this do? it iterates through each cluster
        antigen_counts = Counter([antigen for _, antigen in tcr_list]) #what does this do? it counts the number of times each antigen appears in the cluster
        most_common_antigen, count = antigen_counts.most_common(1)[0] #what does this do? it finds the most common antigen in the cluster
        cluster_size = len(tcr_list) #what does this do? it finds the number of TCRs in the cluster
        sum_most_common += count #what does this do? it adds the number of times the most common antigen appears in the cluster to the sum of most common antigens

        # Determine the purity label
        if count == cluster_size:
            cluster_purity_dict[cluster_id] = "Pure"
        else:
            cluster_purity_dict[cluster_id] = "Impure"

    # Calculate purity and retention rate based on new definitions
    purity = sum_most_common / total_tcrs
    retention_rate = total_tcrs / total_reference_tcrs

    print(f"Purity: {purity}")
    print(f"Retention Rate: {retention_rate}")

    # Output the modified data with the "Cluster Purity" column
    output_file = "clusterPurityResultsV2.tsv"
    with open(clustered_file_path, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Write metadata lines of results of this script
        writer.writerow(["## The number of pure clusters: " + str(sum_most_common)])
        writer.writerow(["## Purity: " + str(purity)])
        writer.writerow(["## Retention Rate: " + str(retention_rate)])

        for row in reader:
            if row[0].startswith('##'):
                writer.writerow(row)  # Write metadata lines from the input file of clustering results
                continue
            cluster_id = row[1]
            purity_label = cluster_purity_dict.get(cluster_id, "Impure")
            writer.writerow(row + [purity_label])

    # Return purity and retention rate
    return purity, retention_rate

# File paths to the clustered output and filtered reference input
clustered_file_path = 'results/filteredTCRantigenData_unique--RotationEncodingBL62.txt'
reference_file_path = 'results/filteredTCRantigenData_unique.txt'
purity, retention_rate = calculate_purity_and_retention(clustered_file_path, reference_file_path)