import csv
import os
import argparse
from collections import Counter

def find_file_in_pwd(filename):
    """Search for a file in the current working directory."""
    for root, dirs, files in os.walk(os.getcwd()):
        if filename in files:
            return os.path.join(root, filename)
    return None

def count_unique_meta_labels(file_path):
    """Count unique meta labels in the specified file."""
    # Initialize a counter for meta labels
    meta_label_counter = Counter()

    # Open the input file and read the data
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) < 6:
                continue  # Skip rows that do not have enough columns
            meta_label = row[5]  # Get the meta label from the 6th column
            meta_label_counter[meta_label] += 1

    # Output the count of unique meta labels
    print(f"Total unique meta labels: {len(meta_label_counter)}")
    print("Meta label counts:", meta_label_counter)

if __name__ == "__main__":
    # Setup argparse to handle command-line arguments
    parser = argparse.ArgumentParser(description="Count unique meta labels in a TSV file.")
    parser.add_argument('-f', '--file', type=str, required=True, help="Name of the input file.")
    args = parser.parse_args()

    # Search for the file in the current working directory
    file_path = find_file_in_pwd(args.file)
    if file_path:
        count_unique_meta_labels(file_path)
    else:
        print(f"File '{args.file}' not found in the current working directory.")