#!/bin/bash                    

#SBATCH --partition=normal    # set the partition where the job will run
#SBATCH --nodes=1             # set the number of nodes 
#SBATCH --cpus-per-task=8     # cores, to run multithreading software
#SBATCH --job-name=CaastTest    # Job name identification
#SBATCH --output=Out.%j.txt        # Path to stdout file
#SBATCH --error=Err.%j.txt         # Path to stderr file
#SBATCH --time=06:59:00  # time (D-HH:MM)
#SBATCH --mail-user=noelia.rodriguez@upf.edu  # Send an email if the job
#SBATCH --mail-type=END,FAIL                   #+  / ends / fails

#Define modules
module load Anaconda3
source activate /homes/users/nrodriguez/.conda/envs/caastoolENV

#Define PATH arguments
BIN=$HOME/caastools/ct
PROTEIN_ALIGNMENTS=$HOME/caastools/examples/MSA/primates.msa.pr
CONFIG_FILE=$HOME/caastools/examples/config.tab
phylogenetic_tree=$HOME/caastools/examples/phylogeny.nw
FAMILIES=$HOME/caastools/test/sp2fam.210727.tab
OUTDIR=$HOME/scratch/caastools_test/out/

# Create directory only when it doesn't exist
if [ ! -d ${OUTDIR} ];then mkdir -p ${OUTDIR};fi

# Run discovery
${BIN} discovery --traitfile ${CONFIG_FILE} \
-a ${PROTEIN_ALIGNMENTS} --fmt phylip-relaxed \
--patterns=1,2 \
-o ${OUTDIR}discovery.caas_test \
--max_fg_gaps=0 \
--max_bg_gaps=0 \
--max_fg_miss=0 \
--max_bg_miss=0

# Run resample
${BIN} resample  \
-p ${phylogenetic_tree}   \
--bytemp ${CONFIG_FILE}  \
-m random --limit_by_group ${FAMILIES} -o ${OUTDIR}phylogeny_restricted.resampling.tab

# Run bootstrap
${BIN} bootstrap  \
-s ${OUTDIR}phylogeny_restricted.resampling.tab  \
-t ${CONFIG_FILE}  \
-a ${PROTEIN_ALIGNMENTS} --fmt phylip-relaxed -o ${OUTDIR}phylogeny_restricted.bootstrap.tab