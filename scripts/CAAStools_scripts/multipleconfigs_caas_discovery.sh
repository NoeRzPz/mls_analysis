#!/bin/bash                   

#SBATCH --array=1-16133%200        # we limit number simultaneously running tasks from the job array to 200
#SBATCH --job-name=Caastool        # Job name identification
#SBATCH --partition=haswell
#SBATCH --output=Out.txt        # Path to stdout file
#SBATCH --error=Err.txt         # Path to stderr file
#SBATCH --time=07:59:00            # time (D-HH:MM)
#SBATCH --mail-user=noelia.rodriguez@upf.edu  # Send an email if the job
#SBATCH --mail-type=END,FAIL                   #+  / ends / fails

#Define modules
module load Anaconda3
source activate /homes/users/nrodriguez/.conda/envs/caastoolENV

#Define PATH arguments
BIN=$HOME/caastools/ct
PROTEIN_ALIGNMENTS=$HOME/scratch/caastools_LQ/data/primate.alignments.231026/
CONFIG_FILES=$HOME/scratch/caastools_LQ/data/2vs2_families/
phylogenetic_tree=$HOME/scratch/caastools_LQ/data/primates.tree.nwk
OUTDIR=$HOME/scratch/caastools_LQ_out/

# Create directory only when it doesn't exist
if [ ! -d ${OUTDIR} ];then mkdir -p ${OUTDIR};fi

# Define arguments in each task
FILE=$(ls ${PROTEIN_ALIGNMENTS}*.phy | awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id')
    
gene_name=$(basename $FILE | cut -d '.' -f 1)

for filepath in $(ls ${CONFIG_FILES}*txt)
do 
    contrast_name=$(basename $filepath | cut -d '.' -f 1)
    if [ ! -d ${OUTDIR}${contrast_name}/ ];then mkdir -p ${OUTDIR}${contrast_name}/;fi

    # Print info of the task
    echo "${SLURM_ARRAY_TASK_ID}, Running CAAStools on ${gene_name} ..." 
    
    # Run discovery  
    ${BIN} discovery --traitfile ${filepath} \
    -a ${FILE} --fmt phylip-relaxed \
    --patterns=1,2,3 \
    -o ${OUTDIR}${contrast_name}/${gene_name}.caas_discovery

done