# Import of Python packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Function to read and process expression data
def read_dataset(batch):
    file_path = Path("Batch") / batch / "expressions.tsv.gz"
    data = pd.read_csv(file_path, sep='\t', compression='gzip')
    return data

# Set up the batches
batches = ["RNAseq1","RNAseq2","BRBseq1","BRBseq2"]

# Perform PCA/gene expression correlation analysis to reveal outliers

def PCA_outlier_removal(data, batch, pca_std_threshold=3):
    # Initialize results dataframe
    results_df = pd.DataFrame(index=data.index)
    results_df['Batch'] = batch
    results_df['Sample_ID'] = results_df.index
    results_df['PCA_Status'] = 'OK'
    # results_df['Correlation_Status'] = 'OK'
    results_df['Final_Status'] = 'OK'
    # Standardize data for PCA
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)
    # Perform PCA
    pca = PCA()
    pca_result = pca.fit_transform(scaled_data)
    # Check PCA outlier
    pc1_mean, pc1_std = np.mean(pca_result[:, 0]), np.std(pca_result[:, 0])
    pc2_mean, pc2_std = np.mean(pca_result[:, 1]), np.std(pca_result[:, 1])
    # Mark PCA outliers
    for idx, (pc1, pc2) in enumerate(pca_result[:, :2]):
        if (abs(pc1 - pc1_mean) > pca_std_threshold * pc1_std or 
            abs(pc2 - pc2_mean) > pca_std_threshold * pc2_std):
            results_df.loc[data.index[idx], 'PCA_Status'] = 'Failed'
    # Set final status
    results_df.loc[(results_df['PCA_Status'] == 'Failed') , 'Final_Status'] = 'Failed'
    # Add PC1 and PC2 values
    results_df['PC1'] = pca_result[:, 0]
    results_df['PC2'] = pca_result[:, 1]
    return results_df

# Process all batches
all_results = []
for batch in batches:
    # Read and transpose data
    data = read_dataset(batch).T
    if data.mean().max() > 20:
        data = np.log2(1+data)
    # Analyze samples
    batch_results = PCA_outlier_removal(data, batch)
    all_results.append(batch_results)
    # Print summary for this batch
    print(f"\nSummary for {batch}:")
    print(f"Total samples: {len(batch_results)}")
    print(f"Passed samples: {sum(batch_results['Final_Status'] == 'OK')}")
    print(f"Failed samples: {sum(batch_results['Final_Status'] == 'Failed')}")

    # PCA plot colored by status
    plt.figure(figsize=(5, 3))
    sns.scatterplot(data=batch_results, x='PC1', y='PC2', hue='Final_Status', style='Final_Status')
    plt.title(f'{batch} - PCA Plot')
    plt.tight_layout()
    plt.savefig(f'./Result/PCA_outliers_detection/{batch}_PCA_outlier_detection.png', dpi=300, bbox_inches='tight')
    plt.close()

    # PCA plot after PCA outliers removal
    plt.figure(figsize=(5, 3))
    sns.scatterplot(data = batch_results[batch_results.Final_Status == "OK"], x='PC1', y='PC2', hue='Final_Status', style='Final_Status')
    plt.title(f'{batch} - PCA Plot (outliers removal)')
    plt.tight_layout()
    plt.savefig(f'./Result/PCA_outliers_removal/{batch}_PCA_outlier_removal.png', dpi=300, bbox_inches='tight')
    plt.close()

# Combine all results
final_results = pd.concat(all_results)
# Save the complete results
final_results.to_csv('sample_qc_results.tsv', sep='\t', index=False)
# Create a summary of results
summary = final_results.groupby(['Batch', 'Final_Status']).size().unstack(fill_value=0)
print("\nOverall Summary:")
print(summary)


