from Bio import AlignIO
import pandas as pd
import os
import re
import glob
import argparse

# Define our data directory
dataDir = '~/mls_analysis/data/'

# Set up the argument parser
parser = argparse.ArgumentParser(description='Process species names.')
parser.add_argument('species', nargs='+', help='List of species names separated by space')

# Parse the command line arguments
args = parser.parse_args()

# Assign the parsed species names to the LQ_spp variable
LQ_spp = args.species
# Define the top species of interest
#LQ_spp = ['Homo_sapiens', 'Macaca_fuscata', 'Macaca_mulatta', 'Macaca_nemestrina', 'Cheirogaleus_medius']

# File path for the dataframe saved from R
file_path = os.path.expanduser("~/mls_analysis/out/functional/df_byFGN.txt")

# Read the discovered CAAS dataframe into a pandas DataFrame
df_byFGN = pd.read_csv(file_path, sep="\t")


# Expand the user's home directory in the dataDir path
expanded_dataDir = os.path.expanduser(dataDir)

# Function to process substitutions with multiple amino acids
def process_substitution(substitution):
    # Extract all letters before the slash using regex
    amino_acids = re.findall(r"([A-Z]+)/", substitution)[0]
    return list(amino_acids) if len(amino_acids) > 1 else [amino_acids]


# This list will hold dictionaries for each finding
results = []

# Check alignment by alignment 
# Only check in multiple alignments in which was found a CAAS in the gene
# Iterate over the genes
for gene in df_byFGN['Gene'].unique():
    # Extract substitutions for the current gene
    gene_data = df_byFGN[df_byFGN['Gene'] == gene]
    # Initialize the substitutions dictionary
    substitutions_of_interest = {}
    for idx, row in gene_data.iterrows():
        # If the position is already in the dictionary, append the amino acid(s)
        position = row['Position'] 
        amino_acids = process_substitution(row['Substitution'])
        
        if position in substitutions_of_interest:
            substitutions_of_interest[position] += amino_acids
        else:
            substitutions_of_interest[position] = amino_acids
    
    # Read the alignment file
    pattern = os.path.join(expanded_dataDir, f"palignments.231026/{gene}.*.phy")
    # Use glob to find files that match the pattern
    # Use glob to find files that match the pattern
    matching_files = glob.glob(pattern)

    # Now you can iterate over all matching files
    for file_path in matching_files:
        try:
            alignment = AlignIO.read(file_path, 'phylip-relaxed')
            # Process the alignment as before
        except Exception as e:
            print(f"An error occurred while processing the file {file_path}: {e}")
    
    
    # Check for each substitution in the species of interest
    for record in alignment:
        if record.id in LQ_spp:
            for position, amino_acids in substitutions_of_interest.items():
                
                # Check if the sequence at the position matches any of the amino acids
                if record.seq[position] in amino_acids:
                    print(f"CAAS {position} {df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0]} in {gene}, is present in {record.id}")
                #else:
                    #print(f"No coincidence for CAAS {position} {df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0]} in {gene}")
                    # Collect the information in a dictionary
                    result = {
                        "Position": position,
                        "Substitution": df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0],
                        "Gene": gene,
                        "Species": record.id
                        }
                    results.append(result)
    
# Convert the list of dictionaries to a DataFrame
results_df = pd.DataFrame(results)

# Define the file path where you want to save the tab-separated file
output_file_path = os.path.expanduser("~/mls_analysis/out/functional/CAAS_validation_results.tsv")

# Save the DataFrame to a tab-separated file
results_df.to_csv(output_file_path, sep='\t', index=False)

print(f"Results have been saved to {output_file_path}")