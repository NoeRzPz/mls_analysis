#!/bin/bash                     #Bash directive

#Define PATH arguments
PROTEIN_ALIGNMENTS=~/mls_analysis/data/example/primates.msa.pr
CONFIG_FILE=~/mls_analysis/data/example/config.tab
PHYLOGENETIC_TREE=~/mls_analysis/data/example/phylogeny.nw
FAMILIES=~/mls_analysis/data/example/sp2fam.210727.tab
OUTDIR=~/mls_analysis/out
mkdir -p ${OUTDIR}

# Run discovery
ct discovery --traitfile ${CONFIG_FILE} \
-a ${PROTEIN_ALIGNMENTS} --fmt phylip-relaxed \
--patterns=1,2 \
-o ${OUTDIR}/caas_test \
--max_fg_gaps=0 \
--max_bg_gaps=0 \
--max_fg_miss=0 \
--max_bg_miss=0

# Run resample
ct resample  \
-p ${PHYLOGENETIC_TREE}   \
--bytemp ${CONFIG_FILE}  \
-m random --limit_by_group ${FAMILIES} -o ${OUTDIR}/phylogeny_restricted.resampling.tab

# Run bootstrap
ct bootstrap  \
-s ${OUTDIR}/phylogeny_restricted.resampling.tab  \
-t ${CONFIG_FILE}  \
-a ${PROTEIN_ALIGNMENTS} --fmt phylip-relaxed -o ${OUTDIR}/phylogeny_restricted.bootstrap.tab