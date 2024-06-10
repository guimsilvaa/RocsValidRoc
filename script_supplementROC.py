import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

# Function to read and preprocess the ROCS output files
def preprocess_rocs_file(file_path):
    df = pd.read_csv(file_path, sep='\t')
    df['Outcome'] = 1 if 'actives' in file_path else 0
    return df[['TanimotoCombo', 'Outcome']]

# Function to extract query number from file name
def extract_query_number(file_name):
    match = re.search(r'queries10_(\d+)', file_name)
    if match:
        return int(match.group(1))
    return None

# Function to calculate ROC and PR metrics
def calculate_metrics(y_true, y_score):
    fpr, tpr, thresholds = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    # Calculate estimated best threshold
    best_threshold_index = np.argmax(tpr - fpr)
    best_threshold = thresholds[best_threshold_index]

    # BEDROC
    alpha = 20  # Adjust alpha as needed
    bedroc = ((1 / alpha) * np.sum(tpr * np.exp(-alpha * fpr))) / ((1 / alpha) * np.sum(np.exp(-alpha * fpr)))

    return fpr, tpr, roc_auc, best_threshold, bedroc

# Function to plot ROC curve
def plot_roc_curve(fpr, tpr, roc_auc, label, color):
    plt.plot(fpr, tpr, label=f'{label} (AUC = {roc_auc:.2f})', color=color)

# Main function
def main():
    # Get all .rpt files in the current directory
    rpt_files = [file for file in os.listdir() if file.endswith('.rpt')]

    # Initialize lists to store metrics
    query_numbers = []
    auc_values = []
    best_thresholds = []
    bedroc_values = []

    # Initialize colors for different curves
    colors = plt.cm.tab10.colors

    # Plot all ROC curves
    plt.figure(figsize=(10, 6))
    for rpt_file in rpt_files:
        if 'actives' in rpt_file:
            active_df = preprocess_rocs_file(rpt_file)
        elif 'decoys' in rpt_file:
            decoy_df = preprocess_rocs_file(rpt_file)
            # Merge active and decoy dataframes
            merged_df = pd.concat([active_df, decoy_df])

            # Calculate ROC and PR metrics
            fpr, tpr, roc_auc, best_threshold, bedroc = calculate_metrics(
                merged_df['Outcome'], merged_df['TanimotoCombo'])

            # Extract query number from file name
            query_number = extract_query_number(rpt_file)
            if query_number is not None:
                query_numbers.append(query_number)
                auc_values.append(roc_auc)
                best_thresholds.append(best_threshold)
                bedroc_values.append(bedroc)

                # Plot ROC curve with query number as label
                plot_roc_curve(fpr, tpr, roc_auc, f'Query {query_number}', colors[query_number - 1])

    # Plot diagonal line representing random performance
    plt.plot([0, 1], [0, 1], linestyle='--', color='lightgrey', label='Random (AUC = 0.50)')

    # Sort legend entries in descending order based on AUC values
    sorted_indices = sorted(range(len(auc_values)), key=lambda i: auc_values[i], reverse=True)
    sorted_query_numbers = [query_numbers[i] for i in sorted_indices]

    # Update legend
    handles, labels = plt.gca().get_legend_handles_labels()
    sorted_handles = [handles[i] for i in sorted_indices]
    sorted_labels = [f'{label} (AUC = {auc_values[i]:.2f})' for i, label in zip(sorted_indices, sorted_query_numbers)]
    # Find the index of "Random" in the original labels list
    random_index = labels.index('Random (AUC = 0.50)')
    # Move "Random" to the end of the labels list
    updated_labels = sorted_labels[:random_index] + sorted_labels[random_index+1:] + ['Random (AUC = 0.50)']
    updated_handles = sorted_handles[:random_index] + sorted_handles[random_index+1:] + [handles[random_index]]
    plt.legend(updated_handles, updated_labels, loc='lower right')

    # Set plot limits
    plt.xlim([0, 1])
    plt.ylim([0, 1])

    # Add labels to x-axis and y-axis
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    # Save the plot
    plt.savefig('ROC_Curve.png')

    # Show plot
    plt.show()

    # Create a DataFrame with the metrics
    metrics_df = pd.DataFrame({
        'Query Number': query_numbers,
        'AUC': auc_values,
        'Best Threshold': best_thresholds,
        'BEDROC': bedroc_values
    })

    # Save the DataFrame to a CSV file
    metrics_df.to_csv('metrics.csv', index=False)

if __name__ == "__main__":
    main()

