# Longevity analysis

## Exploratory Descriptive Analysis of phenotipcic traits

We create eda_r conda environment , to work with traits_eda.Rmd script

```
mamba env create --file envs/eda_r.yaml
```
We dowonload Newick tree (with PHAST estimated branch lengths from 241-way Zoonomia mammalian alignment, update V2 ):

```
wget https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.phast-242.nh
```

## CAAStools

We have followed the instructions on `https://github.com/linudz/caastools`  to intall it.

We create conda environment with dependencies needed by CAAStools

```
mamba env create --file envs/caast_test.yaml
```

We can run localy the example provided with: 

```
bash scripts/local_caas_test.sh
```
Commands for running internal validation:

```
python scripts/internal_validation.py "functional/df_byFGN.txt" "functional/4fam_other_spp.txt"

python scripts/internal_validation.py "functional/cebi_atel_byFGN.txt" "functional/cebi_atel_other_spp.txt"

python scripts/internal_validation.py "functional/lemu_atel_byFGN.txt" "functional/lemu_atel_other_spp.txt"

python scripts/internal_validation.py "functional/cerco_atel_byFGN.txt" "functional/cerco_atel_other_spp.txt"

python scripts/internal_validation.py "functional/cerco_lemu_byFGN.txt" "functional/cerco_lemu_other_spp.txt"

python scripts/internal_validation.py "functional/cerco_cebi_byFGN.txt" "functional/cerco_cebi_other_spp.txt"

python scripts/internal_validation.py "functional/cebi_lemu_byFGN.txt" "functional/cebi_lemu_other_spp.txt"
```
Download mammals protein msa:

```
wget -r -np -nH --cut-dirs=3 --no-check-certificate -A '*.fasta.gz' https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments/

```

Commands for running external validation:

```
python scripts/external_validation.py "df_byFGN.txt" "external_validation/labels.txt"

python scripts/external_validation.py "cebi_atel_byFGN.txt" "external_validation/labels.txt"

python scripts/external_validation.py "cebi_lemu_byFGN.txt" "external_validation/labels.txt"

python scripts/external_validation.py "lemu_atel_byFGN.txt" "external_validation/labels.txt"

python scripts/external_validation.py "cerco_atel_byFGN.txt" "external_validation/labels.txt"

python scripts/external_validation.py "cerco_lemu_byFGN.txt" "external_validation/labels.txt"

python scripts/external_validation.py "cerco_cebi_byFGN.txt" "external_validation/labels.txt"
```
We check the 6 external validated CAAS in our primates, to see if AA substituted are more abundant in primates
```
python scripts/look4external_validated.py
```

We run first Alejandro's file
```
python DATA_RAW_ALIGNMENTS_AND_POSITIONS/Extract_genomic_coordinates.py  DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_cds_files DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.cds.fa.gz DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_html_files DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.sorted.gff out/functional/cebi_lemu_byFGN.txt out/coordinates_cebi_lemu_genes005.tab DATA_RAW_ALIGNMENTS_AND_POSITIONS/proteiID_gene_equivalences.txt out/functional/genes_cebi_lemu005.tsv

python DATA_RAW_ALIGNMENTS_AND_POSITIONS/originalExtract_genomic_coordinates.py  DATA_RAW_ALIGNMENTS_AND_POSITIONS/Complete_CDS_alignments.tar.gz DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.cds.fa.gz DATA_RAW_ALIGNMENTS_AND_POSITIONS/FILTER_CODONS-001.tar.gz DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.sorted.gff out/functional/cebi_lemu_byFGN.txt output_coordinates.tab DATA_RAW_ALIGNMENTS_AND_POSITIONS/proteiID_gene_equivalences.txt ESPN
```
I have downloaded Homo_sapiens.CDS.fa.gz made with the Ensemble 104  annotation version from 
```
https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cds/
```
Info of this file on https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cds/README

```
python DATA_RAW_ALIGNMENTS_AND_POSITIONS/select_gene_files.py  DATA_RAW_ALIGNMENTS_AND_POSITIONS/Complete_CDS_alignments.tar.gz DATA_RAW_ALIGNMENTS_AND_POSITIONS/FILTER_CODONS-001.tar.gz DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_html_files DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_cds_files out/functional/genes_cebi_lemu005.tsv
```
We run transvar script with:
```
bash DATA_RAW_ALIGNMENTS_AND_POSITIONS/run_obtain_protposition.sh 
```
we run last script for the mammals coordenates with:

```
python DATA_RAW_ALIGNMENTS_AND_POSITIONS/second_genomic_coord_aa.py  DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_cds_files DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.cds.fa.gz DATA_RAW_ALIGNMENTS_AND_POSITIONS/selected_html_files DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.sorted.gff out/coordinates_cebi_lemu_genes005.tab DATA_RAW_ALIGNMENTS_AND_POSITIONS/proteiID_gene_equivalences.txt out/functional/genes_cebi_lemu005.tsv DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.GRCh38.cds.all.fa.gz DATA_RAW_ALIGNMENTS_AND_POSITIONS/Homo_sapiens.GRCh38.104.gff3.gz out/mammalcoordinates_cebi_lemu_genes005.tab

```
```
python scripts/external_validation.py "external_validation/mammals_coord.txt" "external_validation/labels.txt"
```