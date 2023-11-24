from Bio import AlignIO
import pandas as pd
import os
import re
import glob
import sys
import gzip

# Define your data directories
dataDir = '~/mls_analysis/data/'
outDir = '~/mls_analysis/out/functional/'

# Expand the user's home directory in the directories path
expanded_dataDir = os.path.expanduser(dataDir)
expanded_outDir = os.path.expanduser(outDir)

# Check if the file path argument is provided
if len(sys.argv) < 3:
    print("Please provide a CAAS file path and species file path as arguments.")
    sys.exit(1)

# Assign the relative file path from the first command-line argument
relative_file_path = sys.argv[1]

# File paths for the dataframes saved from R
# File for families contrast
# Combine the base directory with the relative file path
caasfile_path = os.path.join(expanded_outDir, relative_file_path)

traitfile_path = os.path.join(expanded_outDir, "external_validation/mammals.txt")

# Read the dataframe into a pandas DataFrame
df_byFGN = pd.read_csv(caasfile_path, sep="\t")
traits = pd.read_csv(traitfile_path, sep='\t')

# Function to process substitutions with multiple amino acids
def long_substitution(substitution):
    # Extract all letters before the slash using regex
    long_aa = re.findall(r"([A-Z]+)/", substitution)[0]
    return list(long_aa) if len(long_aa) > 1 else [long_aa]

def short_substitution(substitution):
    # Extract all letters after the slash using regex
    matches = re.findall(r"[A-Z]+/([A-Z]+)", substitution)
    if not matches:
        # Handle the case where there is no match
        return None
    short_aa = matches[0]  # Get the first (and likely only) match
    return list(short_aa) if len(short_aa) > 1 else [short_aa]

# Define the species and the substitutions of interest

# Assign the other species file path from the command-line argument
other_spp_filepath = sys.argv[2]
other_spp_file_path = os.path.join(expanded_outDir,other_spp_filepath)

# Initialize an empty list to store the species names
mammal_spp = []

# Open the file and read each line
with open(other_spp_file_path, 'r') as file:
    for line in file:
        # Add the species name to the list, stripping any trailing newline characters
        mammal_spp.append(line.strip())


# Create a list to hold the data for the new DataFrame
species_classification_data = []


# Check alignment by alignment 
# Only iterating over multiple alignments from genes in which a CAAS was found a CAAS
for gene in df_byFGN['Gene'].unique():
    # Extract substitutions for the current gene
    gene_data = df_byFGN[df_byFGN['Gene'] == gene]
    
    # Initialize the substitutions dictionaries
    substitutions_long = {}
    substitutions_short = {}

    for idx, row in gene_data.iterrows():
        # If the position is already in the dictionary, append the amino acid(s)
        position = row['Position']

        long_aas = long_substitution(row['Substitution'])
        if position in substitutions_long:
            substitutions_long[position] += long_aas
        else:
            substitutions_long[position] = long_aas

        short_aas = short_substitution(row['Substitution'])
        if position in substitutions_short:
            substitutions_short[position] += short_aas
        else:
            substitutions_short[position] = short_aas
            
    # Read the alignment file
    pattern = os.path.join(expanded_outDir, f"external_validation/MultipleCodonAlignments/ENST*.{gene}.fasta.gz")
    
    # Use glob to find files that match the pattern
    matching_files = glob.glob(pattern)

    # Now you can iterate over all matching files
    for file_path in matching_files:
        try:
            # Open the gzip-compressed file in text mode
            with gzip.open(file_path, 'rt') as handle:
                # Read the alignment from the decompressed file
                alignment = AlignIO.read(handle, 'fasta')
        except Exception as e:
            print(f"An error occurred while processing the file {file_path}: {e}")
    
    # Check for each substitution in the species of interest
    for record in alignment:
        if record.id in mammal_spp:
            for position in gene_data['Position'].unique():
                # Default classification
                life_classification = None
                lq = None

                # Check if the adjusted position is within the sequence range
                if position >= 0 and position < len(record.seq):

                    # Check long substitutions
                    if position in substitutions_long and record.seq[position] in substitutions_long[position]:
                        print(f"CAAS {position} {df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0]} in {gene}, {traits.loc[traits['TOGA_ID'] == record.id, 'SPECIES'].values[0]} has long-life AA")
                        life_classification = 'Long'
                        lq = traits.loc[traits['TOGA_ID'] == record.id, 'LQ'].values[0]

                    # Check short substitutions
                    if position in substitutions_short and record.seq[position] in substitutions_short[position]:
                        print(f"CAAS {position} {df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0]} in {gene}, {traits.loc[traits['TOGA_ID'] == record.id, 'SPECIES'].values[0]} has short-life AA")
                        life_classification = 'Short'
                        lq = traits.loc[traits['TOGA_ID'] == record.id, 'LQ'].values[0]

                    # Ensure we only add records where a classification was made
                    if life_classification:
                        # Add the data to our list
                        species_classification_data.append({
                            'Gene': gene,
                            'Position': position,
                            'Substitution': row['Substitution'],
                            'label': traits.loc[traits['TOGA_ID'] == record.id, 'label'].values[0],
                            'type_LQ': life_classification,
                            'LQ': lq
                            })
# Create a DataFrame from the collected data
df_classification = pd.DataFrame(species_classification_data)

# Get the prefix from the first command-line argument
file_name_parts = os.path.basename(sys.argv[1]).split('_')
if len(file_name_parts) >= 2:
    file_prefix = '_'.join(file_name_parts[:2])
else:
    # Handle case where there's only one part or none
    file_prefix = file_name_parts[0] if file_name_parts else "default"

# Append the prefix to the output file name
output_file_name = f"{file_prefix}_sppclassification_results.tab"

# Define the path for the output CSV file
output_file_path = os.path.join(expanded_outDir, f"external_validation/{output_file_name}")

# Export the DataFrame to a tab-separated file
df_classification.to_csv(output_file_path, sep='\t', index=False)

print(f"Results have been saved to {output_file_path}")