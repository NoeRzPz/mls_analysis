#!/bin/bash                   

#SBATCH --array=1-16133%200        # we limit number simultaneously running tasks from the job array to 200
#SBATCH --job-name=Caastool        # Job name identification
#SBATCH --output=Out_%A.txt        # Path to stdout file
#SBATCH --error=Err_%A.txt         # Path to stderr file
#SBATCH --time=06:59:00            # time (D-HH:MM)
#SBATCH --mail-user=noelia.rodriguez@upf.edu  # Send an email if the job
#SBATCH --mail-type=END,FAIL                   #+  / ends / fails

#Define modules
module load Anaconda3
source activate /homes/users/nrodriguez/.conda/envs/caastoolENV

#Define PATH arguments
BIN=$HOME/caastools/ct
PROTEIN_ALIGNMENTS=$HOME/scratch/caastools_LQ/data/primate.alignments.231026/
CONFIG_FILE=$HOME/scratch/caastools_LQ/data/LQ_config.txt
phylogenetic_tree=$HOME/scratch/caastools_LQ/data/primates.tree.nwk
OUTDIR=$HOME/scratch/caastools_LQ_out/

# Create directory only when it doesn't exist
if [ ! -d ${OUTDIR} ];then mkdir ${OUTDIR};fi

# Define arguments in each task
FILE=$(ls ${PROTEIN_ALIGNMENTS}*.phy | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id')
   
gene_name=$(basename $FILE | cut -d '.' -f 1)

# Print info of the task
echo "${SLURM_ARRAY_TASK_ID}, Running CAAStools bootstraping on ${gene_name} ..." >> $HOME/scratch/caastools_LQ/genes_analized.txt
    

# Run bootstrap
${BIN} bootstrap  \
-s ${OUTDIR}phylogeny_restricted.resampling.tab  \
-t ${CONFIG_FILE}  \
-a ${FILE} --fmt phylip-relaxed -o ${OUTDIR}${gene_name}brownian_motion.bootstrap.tab

