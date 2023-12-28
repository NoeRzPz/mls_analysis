# Usage: bash scripts/merge_genes_caas.sh LQ_cebidae_atelidae LQ_cercopi_lemuridae LQ_cebidae_lemuridae LQ_lemuridae_atelidae LQ_cercopi_atelidae LQ_cebi_ateli_vs_lemu_cheiro LQ_cercopi_cebidae

# Store all positional arguments (Directory where we have caas_discovery files we want to merge) in a variable
dirs=$@
#chmod -R u+r ${dir}

for dir in $dirs
do
    echo "Merging files in $dir"
    
    # Take the header from the first file
    head -n 1 $(ls out/caas/${dir}/*.caas_discovery | head -n 1) > out/caas/${dir}/all.caas_discovery.tsv

    # Append the content of all files without headers
    for file in $(ls out/caas/${dir}/*.caas_discovery); 
    do
        tail -n +2 "$file" >> out/caas/${dir}/all.caas_discovery.tsv
    done

done

# Command used in the cluster
#ls scratch/caastools_LQ/data/primate.alignments.231026/ | cut -d '.' -f 1  >> scratch/caastools_LQ/background_genes.txt
echo "Done"
