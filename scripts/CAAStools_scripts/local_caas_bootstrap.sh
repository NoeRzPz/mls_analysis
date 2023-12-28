#!/bin/bash 

#Define PATH arguments
#BIN=ct
PROTEIN_ALIGNMENTS=~/mls_analysis/data/palignments.231026/
RESAMPLIN_FILE=~/mls_analysis/out/caas/caastools_LQ_out_resampling/brownian_motion.resampling.tab
CONFIG_FILE=~/mls_analysis/data/caas_tool_inputs/2vs2_families/LQ_cebidae_lemuridae.txt
OUTDIR=~/mls_analysis/out/caas/

# Create directory only when it doesn't exist
if [ ! -d ${OUTDIR}bootstrap_results/ ];then mkdir ${OUTDIR}bootstrap_results/;fi

# Define arguments in each task
for filepath in $(ls ${PROTEIN_ALIGNMENTS}*.phy)
do 
    gene_name=$(echo $filepath | awk -F'/' '{print $NF}'| cut -d '.' -f 1) 
    
    # Check if the output file for this gene already exists, skip to the next if it does
    if [ -f ${OUTDIR}bootstrap_results/${gene_name}brownian_motion.bootstrap.tab ]; then
        echo "Output for ${gene_name} already exists, skipping..."
        continue
    fi
    # Print info of the task
    echo "Running CAAStools bootstraping on ${gene_name}...">> ${OUTDIR}bootstrap_results/genes_analized_in_bootstrap.txt
    
    # Run bootstrap
    ct bootstrap \
    -s ${RESAMPLIN_FILE}  \
    -t ${CONFIG_FILE}  \
    -a ${filepath} --fmt phylip-relaxed -o ${OUTDIR}bootstrap_results/${gene_name}brownian_motion.bootstrap.tab
    echo
done
#cat ${OUTDIR}bootstrap_results/*brownian_motion.bootstrap.tab > ${OUTDIR}/bootstrap_results/all.boot
#rm ${OUTDIR}bootstrap_results/*brownian_motion.bootstrap.tab
#mv out/caas/*brownian_motion.bootstrap.tab out/caas/bootstrap_results/