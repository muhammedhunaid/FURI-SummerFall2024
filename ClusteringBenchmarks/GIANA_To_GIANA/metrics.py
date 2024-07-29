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
    
    # Calculate purity and retention
    for cluster_id, tcr_list in clusters.items():
        antigen_counts = Counter([antigen for _, antigen in tcr_list])
        most_common_antigen, count = antigen_counts.most_common(1)[0]
        # print(count)
        cluster_size = len(tcr_list)
        purity = count / cluster_size
        
        if purity == 1:
            print("Cluster "+ cluster_id +" is pure")
            pure_clusters += 1
            tcrs_in_pure_clusters += cluster_size
    
    retention_rate = tcrs_in_pure_clusters / total_tcrs
    
    print ("The number of pure clusters:", pure_clusters)
    return pure_clusters / len(clusters), retention_rate

# File path to the clustering output
file_path = 'results/TCRantigenData_unique--RotationEncodingBL62.txt'
purity, retention_rate = calculate_purity_and_retention(file_path)
print(f'Purity(%): {purity*100}')
print(f'Retention Rate: {retention_rate*100}')