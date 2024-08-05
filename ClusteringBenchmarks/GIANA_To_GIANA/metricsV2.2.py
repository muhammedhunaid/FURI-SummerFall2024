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
            antigen = row[-1]
            clusters[cluster_id].append((cdr3_seq, antigen))
            total_tcrs += 1

    sum_most_common = 0
    sum_most_common_alt = 0
    total_tcrs_alt = 0
    cluster_purity_dict = {}
    count_pure_clusters_alt = 0

    # Calculate purity, alternative purity, and retention, and output the results
    for cluster_id, tcr_list in clusters.items():
        antigen_counts = Counter([antigen for _, antigen in tcr_list])
        most_common_antigen, count = antigen_counts.most_common(1)[0]
        cluster_size = len(tcr_list)
        sum_most_common += count

        # Calculate local purity
        local_purity = count / cluster_size

        # # Determine the purity label
        # if local_purity > 0.9:
        #     cluster_purity_dict[cluster_id] = "Pure"
        # else:
        #     cluster_purity_dict[cluster_id] = "Impure"

        # Check for Alternative Purity condition
        if local_purity > 0.9:
            sum_most_common_alt += count
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

    #calculate sensitivity and specificity
    specificity = count_pure_clusters_alt / len(clusters)
    sensitivity = sum_most_common_alt / total_reference_tcrs

    print(f"Purity: {purity}")
    print(f"Retention Rate: {retention_rate}")
    print(f"Purity (Alternative): {purity_alt}")
    print(f"Retention Rate (Alternative): {retention_rate_alt}")
    print(f"Sensitivity: {sensitivity}")
    print(f"Specificity: {specificity}")

    # Output the modified data with the "Cluster Purity" column
    output_file = "hc10s10clusterPurityResultsV2.2.tsv"
    with open(clustered_file_path, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Write metadata lines of results of this script
        writer.writerow(["## The number of most common tcrs in every cluster: " + str(sum_most_common)])
        writer.writerow(["## The number of most common tcrs in every cluster (Alternative): " + str(sum_most_common_alt)])
        writer.writerow(["## The number of total tcrs in every cluster (Alternative): " + str(total_tcrs_alt)])
        writer.writerow(["## Purity: " + str(purity)])
        writer.writerow(["## Retention Rate: " + str(retention_rate)])
        writer.writerow(["## Purity (Alternative): " + str(purity_alt)])
        writer.writerow(["## Retention Rate (Alternative): " + str(retention_rate_alt)])
        writer.writerow(["## Sensitivity: " + str(sensitivity)])
        writer.writerow(["## Specificity: " + str(specificity)])

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
clustered_file_path = 'dataGIANAOriginal/hc10s10--RotationEncodingBL62.txt'
reference_file_path = 'dataGIANAOriginal/hc10s10.txt'
purity, purity_alt, retention_rate, retention_rate_alt = calculate_purity_and_retention(clustered_file_path, reference_file_path)