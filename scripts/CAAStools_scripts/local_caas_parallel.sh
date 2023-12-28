#!/bin/bash                     

#Define PATH arguments
PROTEIN_ALIGNMENTS=~/mls_analysis/data/caas_tool_inputs/primate.test.alignments/
CONFIG_FILE=~/mls_analysis/data/caas_tool_inputs/LQ_config.txt
PHYLOGENETIC_TREE=~/mls_analysis/data/caas_tool_inputs/primates.tree.nwk
FAMILIES=~/mls_analysis/data/caas_tool_inputs/my_sp2fam_210727.tab
OUTDIR=~/mls_analysis/out/caast/
BIN=~/tfm/caastools-main/ct

mkdir -p ${OUTDIR}

# Run discovery
for filepath in $(ls ${PROTEIN_ALIGNMENTS}*.phy)
do 
    gene_name=$(echo $filepath | awk -F'/' '{print $NF}'| cut -d '.' -f 1) 
    
    echo "Running CAAStools discovery on ${gene_name}..."
    
    ${BIN} discovery --traitfile ${CONFIG_FILE} \
    -a ${PROTEIN_ALIGNMENTS}${gene_name}.*.phy --fmt phylip-relaxed \
    --patterns=1,2,3 \
    -o ${OUTDIR}${gene_name}.caas_discovery
    echo "Done"
    echo
done



