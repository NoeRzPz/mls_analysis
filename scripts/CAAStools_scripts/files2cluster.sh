#Usage: bash scripts/files2cluster.sh LQ_config.txt primates.tree.nwk

# Store all positional arguments in a variable
args=$@

# Define PATH argument for data files
#in_path=~/mls_analysis/data/caas_tool_inputs
# Define PATH argument for scripts
in_path=~/mls_analysis/scripts
#to upload data
out_dir=/homes/users/nrodriguez/scratch/caastools_test/data/


# Upload CAAStools input files to cluster:
#preguntarles si pudiera hacerlo con loop
for file in $args
do
    scp ${in_path}/$file nrodriguez@marvin.s.upf.edu:${out_dir}
done
