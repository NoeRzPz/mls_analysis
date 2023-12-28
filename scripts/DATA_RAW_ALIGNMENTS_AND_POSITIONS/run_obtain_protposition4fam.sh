#!/bin/bash

#Define PATH arguments
cds_directory=DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_cds_files
fasta_path=DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.cds.fa.gz
tracking_directory=DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_html_files
gff_file=DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.sorted.gff
caas_file=out/functional/4fam_byFGN.txt
output_coordinates=DATA_RAW_ALIGNMENTS_AND_POSITIONS/4fam_coordinates_transvar.tab
gene_equivalences=DATA_RAW_ALIGNMENTS_AND_POSITIONS/proteiID_gene_equivalences.txt
genes=out/functional/Internal_validation/4fam/internval_genes_4fam.tsv

# Get the list of tracking files, CDS files, and CAAS files
tracking_files=("${tracking_directory}"/*.html)
caas_files=("${caas_directory}"/*.tsv)

# Iterate over each line in the genes file using cat
cat "$genes" | while IFS= read -r gene; do
    # Iterate over each tracking file searching for the HTML pattern file
    for tracking_file in "${tracking_directory}"/"${gene}".*.html; do
        if [ -e "$tracking_file" ]; then
            echo $(basename "${tracking_file}")
            tra="${tracking_file}"
        fi
    done
    # Iterate over each cds file searching for the CDS pattern file
    for cds_file in "${cds_directory}"/"${gene}".*.fasta; do
        if [ -e "$cds_file" ]; then
            echo "$(basename "$cds_file")"
            cds="${cds_file}"
        fi
    done
        
    # Run your Python script with the appropriate parameters
    python DATA_RAW_ALIGNMENTS_AND_POSITIONS/Extract_protein_positions_TRANSVAR.py \
        "${cds}" \
        "${fasta_path}" \
        "${tra}" \
        "${gff_file}" \
        "${caas_file}" \
        "${output_coordinates}" \
        "${gene_equivalences}" \
        "${gene}"
done
