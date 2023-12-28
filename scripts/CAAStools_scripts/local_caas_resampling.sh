#!/bin/bash                   

#Define PATH arguments
#BIN=ct
CONFIG_FILE=~/mls_analysis/data/caas_tool_inputs/2vs2_families/LQ_cebidae_lemuridae.txt
phylogenetic_tree=~/mls_analysis/data/caas_tool_inputs/primates.tree.nwk
FAMILIES=~/mls_analysis/data/caas_tool_inputs/my_sp2fam_210727.tab
TRAIT_FILE=~/mls_analysis/data/caas_tool_inputs/LQ.tab
OUTDIR=~/mls_analysis/out/caas/caastools_LQ_out_resampling/

# Create directory only when it doesn't exist
if [ ! -d ${OUTDIR} ];then mkdir ${OUTDIR};fi


# Run resample

# Restricted by family
#${BIN} resample  \
#-p ${phylogenetic_tree}   \
#--bytemp ${CONFIG_FILE}  \
#-m random --limit_by_group ${FAMILIES} -o ${OUTDIR}phylogeny_restricted.resampling.tab


# Brownian motion
ct resample  \
-p ${phylogenetic_tree}   \
--bytemp ${CONFIG_FILE}  \
-m bm --traitvalues ${TRAIT_FILE} -o ${OUTDIR}brownian_motion.resampling.tab

echo "Resampling Done"