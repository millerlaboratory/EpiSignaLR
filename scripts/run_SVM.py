import pickle
import os
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from sklearn import svm
from pathlib import Path

# Adapted from Geysens et al., 2025 from the Vermeesch Lab

# Argument parsing
parser = argparse.ArgumentParser(description="Run SVM classifier on a single BED file sample.")
parser.add_argument('--bed', type=str, required=True, help='Path to the single BED file')
parser.add_argument('--out_prefix', type=str, default='svm_output', help='Prefix for output files')
args = parser.parse_args()

# path to pickle file containing a dictionary with illumina data preprocessed
script_dir = Path(__file__).parent
repo_root = script_dir.parent
file_path = repo_root / 'data' / 'no_strand_all_points_dict.pickle'

with open(file_path, 'rb') as file:
    epi_signatures_all = pickle.load(file)

# Path to the single BED file
bed_file_path = args.bed

# Read the BED file (assuming no header)
bed_df = pd.read_csv(bed_file_path, sep='\t', header=None)
bed_df = bed_df[[0, 1, 2, 3, 14]].reset_index(drop=True)
bed_df = bed_df.rename(columns={0: 'Chr', 1: 'Start', 2: 'End', 3: 'Disorder', 14: 'Methylation'})

# Clean methylation column
bed_df['Methylation'] = bed_df['Methylation'].replace('.', np.nan)
bed_df['Methylation'] = pd.to_numeric(bed_df['Methylation'], errors='coerce')
bed_df['Methylation'] = bed_df['Methylation'] / 100

# Split rows with multiple disorders into separate rows
bed_df = bed_df.assign(Disorder=bed_df['Disorder'].str.split(','))
bed_df = bed_df.explode('Disorder')
bed_df['Disorder'] = bed_df['Disorder'].str.strip()

# Create a dictionary of DataFrames, one per disorder
disorder_dfs = {disorder: df.drop(columns='Disorder').reset_index(drop=True)
                for disorder, df in bed_df.groupby('Disorder')}

# Prepare the reference data with the sample methylation values
epi_signatures_all_with_sample = epi_signatures_all.copy()
epi_signatures_all_with_sample['MRXCJS'] = epi_signatures_all_with_sample['MRXSCJ']
del epi_signatures_all_with_sample['MRXSCJ']

sample_name = os.path.basename(bed_file_path).replace('.bed', '')
print("Matching disorders from sample to reference disorder signatures:")
for disorder in epi_signatures_all_with_sample:
    epi_signatures_all_with_sample[disorder] = epi_signatures_all_with_sample[disorder].iloc[:, 2:-1]
    if disorder in disorder_dfs:
        sample_meth = disorder_dfs[disorder]['Methylation']
        epi_signatures_all_with_sample[disorder][sample_name] = sample_meth.values
        print(f"  - Sample disorder '{disorder}' matched and methylation values added.")
    else:
        print(f"  - Sample disorder '{disorder}' not present in sample BED file; skipping.")

n_of_samples = 1
samples = [sample_name]
print(sample_name)

filter_df = pd.DataFrame()
pred_df = pd.DataFrame()
order_disorders = []

for disorder in epi_signatures_all_with_sample:
    if disorder != 'Controls':
        X = epi_signatures_all_with_sample[disorder].T
        y = X.index.tolist()

        if disorder != 'MRXCJS':
            new_y = [1 if item == str(disorder) else 0 for item in y]
        else:
            new_y = [1 if item == 'MRXSCJ' else 0 for item in y]

        X = X.replace('NA', np.nan)
        X = X.astype(float)
        X = X.fillna(X.mean())

        X_train = X.iloc[:-n_of_samples, :]
        X_test = X.iloc[-n_of_samples:, :]

        y_train = new_y[:-n_of_samples]
        y_test = new_y[-n_of_samples:]

        class_weights = {0: 1, 1: 10}
        linear_classifier = svm.SVC(kernel='linear', class_weight=class_weights)
        linear_classifier.fit(X_train, y_train)

        linear_pred = linear_classifier.predict(X_test)
        new_row_df = pd.DataFrame([linear_pred])
        pred_df = pd.concat([pred_df, new_row_df], ignore_index=True)
        order_disorders.append(disorder)

        decision_values = linear_classifier.decision_function(X_test)
        new_row_df = pd.DataFrame([decision_values])
        filter_df = pd.concat([filter_df, new_row_df], ignore_index=True)

filter_df.index = order_disorders
filter_df.columns = samples

# Output all scores to CSV
scores_csv = f"{args.out_prefix}_scores.csv"
filter_df.to_csv(scores_csv)

# Find top result and write to text file
for sample in range(n_of_samples):
    assigned_values = []
    assigned_disorders = []
    i = 0
    max_value = filter_df.iloc[:, sample].max()
    row_name = filter_df.iloc[:, sample].idxmax()

    for svm_value in filter_df.iloc[:, sample]:
        if svm_value >= 0.35:
            assigned_values.append(svm_value)
            assigned_disorders.append(order_disorders[i])
        i += 1

    # Write top result to text file
    top_txt = f"{args.out_prefix}_top_result.txt"
    with open(top_txt, 'w') as f:
        if max_value == -1000 or max_value < 0.35:
            row_name = 'Control'
        f.write(f"Sample: {samples[sample]}\n")
        f.write(f"Max value: {max_value}\n")
        f.write(f"Assigned Class: {row_name}\n")
        if len(assigned_values) > 1:
            f.write("More than one disorder is detected:\n")
            f.write(str(assigned_disorders) + "\n")
            f.write(str(assigned_values) + "\n")

    # Barplot of top 5 scores
    top5 = filter_df.iloc[:, sample].sort_values(ascending=False).head(5)
    plt.figure(figsize=(8, 5))
    bars = plt.bar(top5.index, top5.values, color='skyblue')
    plt.ylabel('SVM Score')
    plt.title(f'Top 5 SVM Scores for {samples[sample]}')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.2f}', va='bottom', ha='center')
    plt.savefig(f"{args.out_prefix}_top5_barplot.png")
    plt.close()