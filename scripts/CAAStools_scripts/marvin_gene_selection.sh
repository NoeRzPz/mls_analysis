#!/bin/bash                     
#SBATCH --array=1-2191%200      #we limit number simultaneously running tasks from the job array to 100
#SBATCH --job-name=Genesselected   # Job name identification
#SBATCH --output=Out.txt        # Path to stdout file
#SBATCH --error=Err.txt         # Path to stderr file
#SBATCH --time=04:59:00
#SBATCH --mail-user=noelia.rodriguez@upf.edu  # Send an email if the job
#SBATCH --mail-type=END,FAIL                   #+  / ends / fails

#Define modules
module load Anaconda3
source activate /homes/users/nrodriguez/.conda/envs/rercon_caas

#Define PATH arguments
tracking_directory=$HOME/scratch/DATA_RAW_ALIGNMENTS_AND_POSITIONS/FILTER_CODONS-001
cds_directory=$HOME/scratch/DATA_RAW_ALIGNMENTS_AND_POSITIONS/Complete_CDS_alignments
genes=$HOME/scratch/DATA_RAW_ALIGNMENTS_AND_POSITIONS/genes_cebi_lemu005.tsv
html_dir=$HOME/scratch/DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_html_files
cds_dir=$HOME/scratch/DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_cds_files

# Create the directories if they don't exist
if [ ! -d ${html_dir} ];then mkdir ${html_dir};fi
if [ ! -d ${cds_dir} ];then mkdir ${cds_dir};fi

# Define arguments in each task
gene_name=$(awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${genes})

# Print info of the task
echo "${SLURM_ARRAY_TASK_ID}, ${gene_name}"

firstchar=${gene_name:0:1} # First character
secondchar=${gene_name:0:2} # First two characters
  
# Define patterns for HTML and CDS files
html_pattern="${gene_name}.Homo_sapiens.html"
cds_pattern="${gene_name}.Homo_sapiens.fasta"

# Construct the paths to the specific subfolders in tracking_directory and cds_directory
tracking_subfolder="${tracking_directory}/${firstchar}/${secondchar}"
cds_subfolder="${cds_directory}/${firstchar}/${secondchar}"

# Extract HTML files
if [ -d "${tracking_subfolder}" ]; then
    for file in "${tracking_subfolder}"/*.html; do
        if [[ $(basename "$file") == "${html_pattern}" ]]; then
            cp "${file}" "${html_dir}/"
            break
        fi
    done
fi

# Extract CDS files
if [ -d "${cds_subfolder}" ]; then
    for file in "${cds_subfolder}"/*.fasta; do
        if [[ $(basename "$file") == "${cds_pattern}" ]]; then
            cp "${file}" "${cds_dir}/"
            break
        fi
    done
fi
