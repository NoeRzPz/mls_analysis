#Usage: bash scripts/cluster2local.sh Err.21942398.txt Out.21942398.txt
# Store all positional arguments in a variable
args=$@

# Define PATH argument for data files
cluster_path=/homes/users/nrodriguez/scratch

# Local path to download data
local_path=~/mls_analysis/out/caas/

# Upload CAAStools input files to cluster:
#preguntarles si pudiera hacerlo con loop
for file in $args
do
    scp nrodriguez@marvin.s.upf.edu:${cluster_path}/$file ${local_path}
done
