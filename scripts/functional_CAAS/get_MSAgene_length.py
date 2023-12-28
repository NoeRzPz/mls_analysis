from Bio import AlignIO
import pandas as pd
import os
import glob

# Define your data directories
dataDir = '~/mls_analysis/data/'
outDir = '~/mls_analysis/out/'

# Expand the user's home directory in the directories path
expanded_dataDir = os.path.expanduser(dataDir)
expanded_outDir = os.path.expanduser(outDir)


# Create a list to hold the data for the new DataFrame
gene_length_data = []

# Create a set to track which genes have been processed
processed_genes = set()

# Loop through all .phy files in the specified directory
pattern = os.path.join(expanded_dataDir, "palignments.231026/*.phy")

# Use glob to find files that match the pattern
matching_files = glob.glob(pattern)

for file_path in matching_files:
    try:
        # Extract the gene name from the file path
        gene = os.path.basename(file_path).split('.')[0]

        # Read the alignment file
        alignment = AlignIO.read(file_path, 'phylip-relaxed')

        # Process the alignment to get gene length
        for record in alignment:
            gene_length_data.append({
                'Gene': gene,
                'Gene_length_msa': len(record.seq)
            })

    except Exception as e:
        print(f"An error occurred while processing the file {file_path}: {e}")

# Create a DataFrame from the collected data
df_genes_length = pd.DataFrame(gene_length_data)

# Remove rows with duplicated gene
df_genes_length = df_genes_length.drop_duplicates(subset=['Gene'])

# Define the filename for the output CSV file
output_filename = "gene_lengths.tsv"  # You can change this filename as needed

# Define the path for the output CSV file
output_file_path = os.path.join(expanded_outDir, "functional", output_filename)

# Export the DataFrame to a tab-separated file
df_genes_length.to_csv(output_file_path, sep='\t', index=False)

print(f"Results have been saved to {output_file_path}")