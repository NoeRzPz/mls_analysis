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

# File paths for the dataframes saved from R
caasfile_path = os.path.join(expanded_outDir, "functional/df_byFGN.txt")
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

# Define the species and the substitutions of interest
intermediate_quantile_spp = ['Callithrix_jacchus','Cebus_albifrons','Theropithecus_gelada','Mandrillus_sphinx','Mandrillus_leucophaeus','Cercopithecus_neglectus','Cercopithecus_mona','Allenopithecus_nigroviridis','Erythrocebus_patas','Rhinopithecus_roxellana','Colobus_angolensis','Pan_paniscus','Pan_troglodytes','Gorilla_gorilla','Carlito_syrichta','Daubentonia_madagascariensis','Eulemur_fulvus','Eulemur_macaco','Lemur_catta','Microcebus_murinus','Nycticebus_coucang','Otolemur_garnettii']

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
        short_aas = short_substitution(row['Substitution'])

        if (position in substitutions_long) and (position in substitutions_short):
            substitutions_long[position] += long_aas
            substitutions_short[position] += short_aas
        else:
            substitutions_long[position] = long_aas
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
        if record.id in intermediate_quantile_spp:
            # Default classification
            life_classification = None
            lq = None
            for position, amino_acids in substitutions_long.items():
                # Check if the sequence at the position matches any of the amino acids
                if record.seq[position] in amino_acids:
                    print(f"CAAS {position} {df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0]} in {gene}, {record.id} has long-life AA")
                    life_classification = 'Long'
                    lq = traits.loc[traits['label'] == record.id, 'LQ'].values[0]

            for position, amino_acids in substitutions_short.items():
                # Check if the sequence at the position matches any of the amino acids
                if record.seq[position] in amino_acids:
                    print(f"CAAS {position} {df_byFGN.loc[df_byFGN['Position'] == position, 'Substitution'].values[0]} in {gene}, {record.id} has short-life AA")
                    life_classification = 'Short'
                    lq = traits.loc[traits['label'] == record.id, 'LQ'].values[0]

            # Ensure we only add records where a classification was made
            if life_classification:
                # Add the data to our list
                species_classification_data.append({
                    'Gene': gene,
                    'Position': position,
                    'Substitution': row['Substitution'],
                    'label': record.id,
                    'type_LQ': life_classification,
                    'LQ': lq
                    })
# Create a DataFrame from the collected data
df_classification = pd.DataFrame(species_classification_data)

# Define the path for the output CSV file
output_file_path = os.path.join(expanded_outDir, "functional/sppclassification_results.tab")

# Export the DataFrame to a tab-separated file
df_classification.to_csv(output_file_path, sep='\t', index=False)

print(f"Results have been saved to {output_file_path}")