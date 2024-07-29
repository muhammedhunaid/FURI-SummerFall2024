import csv
from collections import defaultdict, Counter

def calculate_purity_and_retention(file_path):
    clusters = defaultdict(list)
    total_tcrs = 0
    
    # Read the file and parse clusters
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('##'):
                continue  # Skip metadata lines
            cdr3_seq = row[0]
            cluster_id = row[1]
            antigen = row[-1].split(':')[-1]
            clusters[cluster_id].append((cdr3_seq, antigen))
            total_tcrs += 1
    
    pure_clusters = 0
    tcrs_in_pure_clusters = 0
    cluster_purity_dict = {}

    # Calculate purity and retention and output the results
    for cluster_id, tcr_list in clusters.items():
        antigen_counts = Counter([antigen for _, antigen in tcr_list])
        most_common_antigen, count = antigen_counts.most_common(1)[0]
        cluster_size = len(tcr_list)
        purity = count / cluster_size
        
        if purity == 1:
            pure_clusters += 1
            tcrs_in_pure_clusters += cluster_size
            cluster_purity_dict[cluster_id] = "Pure"
        else:
            cluster_purity_dict[cluster_id] = "Impure"
    
    retention_rate = tcrs_in_pure_clusters / total_tcrs
    
    print("The number of pure clusters:", pure_clusters)
    print(f'Purity(%): {pure_clusters / len(clusters) * 100}')
    print(f'Retention Rate(%): {retention_rate * 100}')

    # Output the modified data with the "Cluster Purity" column
    output_file = "clusterPurityResults.tsv"
    with open(file_path, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        # Write metadata lines of results of this script
        writer.writerow(["## The number of pure clusters: " + str(pure_clusters)])
        writer.writerow(["## Purity(%): " + str(pure_clusters / len(clusters) * 100)])
        writer.writerow(["## Retention Rate(%): " + str(retention_rate * 100)])

        for row in reader:
            if row[0].startswith('##'):
                writer.writerow(row)  # Write metadata lines from the input file of clustering results
                continue
            cluster_id = row[1]
            purity_label = cluster_purity_dict.get(cluster_id, "Impure")
            writer.writerow(row + [purity_label])

    # Return purity and retention rate
    return pure_clusters / len(clusters), retention_rate

# File path to the clustering output
file_path = 'results/TCRantigenData_unique--RotationEncodingBL62.txt'
purity, retention_rate = calculate_purity_and_retention(file_path)