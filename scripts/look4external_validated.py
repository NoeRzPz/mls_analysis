from Bio import AlignIO
import pandas as pd
import os
import re
import glob

# Define your data directories
dataDir = '~/mls_analysis/data/'
outDir = '~/mls_analysis/out/'

# Expand the user's home directory in the directories path
expanded_dataDir = os.path.expanduser(dataDir)
expanded_outDir = os.path.expanduser(outDir)

# File paths for the external validated CAAS dataframe
# Combine the base directory with the relative file path
caasfile_path = os.path.join(expanded_outDir, "functional/external_validation/external_validated_df.tab")

traitfile_path = os.path.join(expanded_outDir, "primates_traits.csv")

# Read the dataframe into a pandas DataFrame
df_byFGN = pd.read_csv(caasfile_path, sep="\t")
traits = pd.read_csv(traitfile_path)

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

# Define the substitutions of interest
# Create a list to hold the data for the new DataFrame
externalvalidated_caas_primatesdata = []


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
    pattern = os.path.join(expanded_dataDir, f"palignments.231026/{gene}.*.phy")
    
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
        #if record.id in intermediate_quantile_spp:
            for position in gene_data['Position'].unique():
                # Default classification
                life_classification = None
                lq = None

                # Check long substitutions
                if position in substitutions_long and record.seq[position] in substitutions_long[position]:
                    print(f"CAAS {position} {df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0]} in {gene}, {record.id} has long-life AA")
                    life_classification = 'Long'
                    lq = traits.loc[traits['label'] == record.id, 'LQ'].values[0]

                # Check short substitutions
                if position in substitutions_short and record.seq[position] in substitutions_short[position]:
                    print(f"CAAS {position} {df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0]} in {gene}, {record.id} has short-life AA")
                    life_classification = 'Short'
                    lq = traits.loc[traits['label'] == record.id, 'LQ'].values[0]

                # Ensure we only add records where a classification was made
                if life_classification:
                    # Add the data to our list
                    externalvalidated_caas_primatesdata.append({
                        'Gene': gene,
                        'Position': position,
                        'Substitution': row['Substitution'],
                        'label': record.id,
                        'type_LQ': life_classification,
                        'aa': record.seq[position],
                        'LQ': lq
                        })
                else:
                    externalvalidated_caas_primatesdata.append({
                        'Gene': gene,
                        'Position': position,
                        'Substitution': row['Substitution'],
                        'label': record.id,
                        'type_LQ': none,
                        'aa': record.seq[position],
                        'LQ': lq
                        })
# Create a DataFrame from the collected data
df_classification = pd.DataFrame(externalvalidated_caas_primatesdata)

# Append the prefix to the output file name
output_file_name = f"external_validated_CAAS_in_primates.tab"

# Define the path for the output CSV file
output_file_path = os.path.join(expanded_outDir, f"functional/external_validation/{output_file_name}")

# Export the DataFrame to a tab-separated file
df_classification.to_csv(output_file_path, sep='\t', index=False)

print(f"Results have been saved to {output_file_path}")