# Usage: bash scripts/local_transfer_test.sh LQ_config.txt primates.tree.nwk
# Store all positional arguments in a variable
args=$@

#Define PATH argument
in_path=~/mls_analysis/data/caas_tool_inputs
out_dir=~/mls_analysis/

# Upload CAAStools input files to cluster:
#preguntarles si pudiera hacerlo con loop
for file in $args
do
    cp ${in_path}/$file ${out_dir}
done
