# Import of Python packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# Function to read and process expression data
def read_dataset(batch):
    file_path = Path("Batch") / batch / "expressions.tsv.gz"
    data = pd.read_csv(file_path, sep='\t', compression='gzip')
    return data

# Set up the batches
batches = ["RNAseq1","RNAseq2","BRBseq1","BRBseq2"]

for batch in batches:
    # Read data
    data = read_dataset(batch).T
    if data.mean().max() > 20:
        data = np.log2(1+data)
    #Generating distribution plot
    fig, axes = plt.subplots(1, 1, figsize=(5, 3))
    sns.distplot(data.mean(),kde = True)
    axes.grid(False)
    plt.title(f'{batch} - Data distribution check')             
    plt.ylabel('Density')
    plt.savefig(f'./Result/Data_distribution_check/{batch}_Data_distribution_check.png', dpi=300, bbox_inches='tight')
    plt.close()