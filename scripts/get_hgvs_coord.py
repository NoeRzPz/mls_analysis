import requests
import os
import pandas as pd

# Define your data directories
outDir = '~/mls_analysis/out/functional/'

# Expand the user's home directory in the directories path
expanded_outDir = os.path.expanduser(outDir)

# Combine the base directory with the relative file path
transvar_path = os.path.join(expanded_outDir, "Internal_validation/cebi_lemu/internal_val_transvar_coord.tsv")

# Read the dataframe into a pandas DataFrame
df_transvar = pd.read_csv(transvar_path, sep="\t")

# List of HGVS coordinates you want to search for
hgvs_coordinates = df_transvar['transvar_coord'].tolist()

# Initialize an empty DataFrame to store results
results_df = pd.DataFrame(columns=['HGVS', 'Position', 'PosName', 'Description'])

# This function will construct the URL for the UCSC API call
def construct_query_url(hgvs_coord):
    # Replace with the actual UCSC API endpoint for HGVS coordinate queries
    base_url = "https://api.genome.ucsc.edu/search?"
    params = {
        'search': hgvs_coord, 
        'genome': 'hg38',  # Adjust if you are querying a different genome version
    }
    # Constructs the full URL with parameters for the request
    query_url = f"{base_url}{'&'.join(f'{key}={value}' for key, value in params.items())}"
    return query_url

# This function will perform the actual API call to UCSC and save results
def query_ucsc(hgvs_coord):
    query_url = construct_query_url(hgvs_coord)
    response = requests.get(query_url)
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()
        # Extract the 'positionMatches' field which is a list
        matches = data.get('positionMatches', [])
        # Iterate through each match and collect relevant data
        for match in matches:
            for m in match.get('matches', []):
                # Append each result as a new row in the DataFrame
                results_df.loc[len(results_df)] = [hgvs_coord, m.get('position'), m.get('posName'), match.get('description')]
        return data
    else:
        # Handle the error here
        print(f"Failed to retrieve data for {hgvs_coord}: {response.status_code}")
        return None

# Iterate over the list of HGVS coordinates and perform the queries
for hgvs_coord in hgvs_coordinates:
    query_ucsc(hgvs_coord)

# Define the path for the output file
output_path = os.path.join(expanded_outDir, "Internal_validation/cebi_lemu/ucsc_results.tsv")

# Save the DataFrame to a TSV file
results_df.to_csv(output_path, sep='\t', index=False)
print(f"Results have been saved to {output_path}")