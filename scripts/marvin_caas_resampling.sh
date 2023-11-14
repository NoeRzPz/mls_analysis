#!/bin/bash                   

#SBATCH --partition=normal    # set the partition where the job will run
#SBATCH --nodes=1             # set the number of nodes 
#SBATCH --cpus-per-task=8     # cores, to run multithreading software
#SBATCH --job-name=Caasresampling        # Job name identification
#SBATCH --output=Out.%j.txt        # Path to stdout file
#SBATCH --error=Err.%j.txt         # Path to stderr file
#SBATCH --time=01:59:00            # time (D-HH:MM)
#SBATCH --mail-user=noelia.rodriguez@upf.edu  # Send an email if the job
#SBATCH --mail-type=END,FAIL                   #+  / ends / fails

#Define modules
module load Anaconda3
source activate /homes/users/nrodriguez/.conda/envs/rercon_caas

#Define PATH arguments
BIN=$HOME/caastools/ct
CONFIG_FILE=$HOME/scratch/caastools_LQ/data/LQ_config.txt
phylogenetic_tree=$HOME/scratch/caastools_LQ/data/primates.tree.nwk
FAMILIES=$HOME/scratch/caastools_LQ/data/my_sp2fam_210727.tab
TRAIT_FILE=$HOME/scratch/caastools_LQ/data/LQ.tab
OUTDIR=$HOME/scratch/caastools_LQ_out/

# Create directory only when it doesn't exist
if [ ! -d ${OUTDIR} ];then mkdir ${OUTDIR};fi


# Run resample

# Restricted by family
#${BIN} resample  \
#-p ${phylogenetic_tree}   \
#--bytemp ${CONFIG_FILE}  \
#-m random --limit_by_group ${FAMILIES} -o ${OUTDIR}phylogeny_restricted.resampling.tab


# Brownian motion
${BIN} resample  \
-p ${phylogenetic_tree}   \
--bytemp ${CONFIG_FILE}  \
-m bm --traitvalues ${TRAIT_FILE} -o ${OUTDIR}brownian_motion.resampling.tab

echo "Resampling Done"
