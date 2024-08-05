import csv
from collections import defaultdict, Counter

#Added feature to count number of TCR attached to multiple most common antigen in each cluster
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
    sum_most_common_alt = 0
    total_tcrs_alt = 0
    cluster_purity_dict = {}
    count_pure_clusters_alt = 0

    # Calculate purity, alternative purity, retention, alternative retention and output the results
    for cluster_id, tcr_list in clusters.items():
        antigen_counts = Counter([antigen for _, antigen in tcr_list])
        most_common_count = antigen_counts.most_common(1)[0][1]

        # Calculate the total count for the most common antigens
        most_common_total = sum(count for _, count in antigen_counts.items() if count == most_common_count)
        print("Most common total for cluster", cluster_id, most_common_total)
        cluster_size = len(tcr_list)
        sum_most_common += most_common_total

        # Calculate local purity
        local_purity = most_common_total / cluster_size

        # Check for Alternative Purity condition
        if local_purity == 1.0:
            sum_most_common_alt += most_common_total
            total_tcrs_alt += cluster_size
            cluster_purity_dict[cluster_id] = "Pure"
            count_pure_clusters_alt += 1
        else:
            cluster_purity_dict[cluster_id] = "Impure"

    # Calculate purity and retention rate based on new definitions
    purity = sum_most_common / total_tcrs
    purity_alt = sum_most_common_alt / total_tcrs_alt if total_tcrs_alt != 0 else 0
    retention_rate = total_tcrs / total_reference_tcrs
    retention_rate_alt = total_tcrs_alt / total_reference_tcrs

    # Calculate sensitivity and specificity
    specificity = count_pure_clusters_alt / len(clusters)
    #sensitivity = total_tcrs_alt / total_reference_tcrs

    print(f"The number of clusters: {str(len(clusters))}")
    print(f"The number of pure clusters (Local Purity == 1): {count_pure_clusters_alt}")
    print(f"Purity (Using TCR's from all clusters): {purity}")
    #print(f"Retention Rate: {retention_rate}")
    print(f"Purity (Using TCR's from Clusters with Local Purity == 1): {purity_alt}")
    print(f"Retention Rate (Using TCR's from Clusters with Local Purity == 1): {retention_rate_alt}")
    print(f"Sensitivity (Retention): {retention_rate_alt}")
    print(f"Specificity (Fraction): {specificity}")

    # Output the modified data with the "Cluster Purity" column
    output_file = "grayClusterPurityResultsV2.1.1.tsv"
    with open(clustered_file_path, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Write metadata lines of results of this script
        writer.writerow(["## The number of clusters: " + str(len(clusters))])
        writer.writerow(["## The number of pure clusters (Local Purity == 1):  " + str(count_pure_clusters_alt)])
        writer.writerow(["## The number of most common tcrs in every cluster: " + str(sum_most_common)])
        writer.writerow(["## The number of most common tcrs in pure clusters (Local Purity == 1): " + str(sum_most_common_alt)])
        writer.writerow(["## The number of total tcrs in pure cluster (Local Purity == 1): " + str(total_tcrs_alt)])
        writer.writerow(["## The number of total tcrs in every cluster: " + str(total_tcrs)])
        writer.writerow(["## The number of total tcrs in every cluster: " + str(total_reference_tcrs)])
        writer.writerow(["## Purity (Using TCR's from all clusters):: " + str(purity)])
        #writer.writerow(["## Retention Rate: " + str(retention_rate)])
        writer.writerow(["## Purity (Using TCR's from Clusters with Local Purity == 1): " + str(purity_alt)])
        writer.writerow(["## Retention Rate (Using TCR's from Clusters with Local Purity == 1): " + str(retention_rate_alt)])
        writer.writerow(["## Sensitivity (Retention): " + str(retention_rate_alt)])
        writer.writerow(["## Specificity (Fraction): " + str(specificity)])

        for row in reader:
            if row[0].startswith('##'):
                writer.writerow(row)  # Write metadata lines from the input file of clustering results
                continue
            cluster_id = row[1]
            purity_label = cluster_purity_dict.get(cluster_id, "Impure")
            writer.writerow(row + [purity_label])

    # Return purity, alternative purity, and retention rate
    return purity, purity_alt, retention_rate, retention_rate_alt

# File paths to the clustered output and filtered reference input
clustered_file_path = 'results/filteredTCRantigenData_unique--RotationEncodingBL62.txt'
reference_file_path = 'results/filteredTCRantigenData_unique.txt'
purity, purity_alt, retention_rate, retention_rate_alt = calculate_purity_and_retention(clustered_file_path, reference_file_path)