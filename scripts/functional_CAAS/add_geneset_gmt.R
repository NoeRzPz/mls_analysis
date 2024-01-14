# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(GSEABase)
library(readxl)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data/aging")
resultsDir <- file.path(workingDir,"out") 

# Load data
GenAge_human <- read_csv(file.path(dataDir,"/genage_human.csv"))
gnomad <- read_delim(file.path(dataDir,"/gnomad.v4.0.constraint_metrics.tsv"), 
           delim = "\t", escape_double = FALSE, trim_ws = TRUE)
GTEx <- read_excel(file.path(dataDir,"/GTEx41420_2018_93_MOESM2_ESM.xlsx"))

# Read GMT file
ageing_gsets <- getGmt(file.path(workingDir,"data/agein.Hs.symbols.gmt"))

# check number of  gene sets this object has. Each gene set contains the smbols for the genes in each gene set
length(ageing_gsets)

# Function to add a gene set to an existing GeneSetCollection
add_gene_set <- function(gene_sets, new_set_name, genes) {
  # Create a new GeneSet object
  new_gene_set <- GeneSet(genes, setName = new_set_name)
  
  # Add the new GeneSet to the GeneSetCollection
  if (is(gene_sets, "GeneSetCollection")) {
    gene_sets <- GeneSetCollection(c(gene_sets, list(new_gene_set)), geneIdType = geneIdType(gene_sets))
  } else {
    # If gene_sets is not a GeneSetCollection, create a new one
    gene_sets <- GeneSetCollection(list(new_gene_set), geneIdType = "EntrezIdentifier")
  }
  
  return(gene_sets)
}

# Prepare data for New gene set data 
new_set_name <- "Human GenAge"
new_genes <- GenAge_human$symbol

# Add new gene set to the list
ageing_gsets <- add_gene_set(ageing_gsets, new_set_name, new_genes)
length(ageing_gsets)
# Verify if the gene set has been added
print(ageing_gsets[[new_set_name]])

new_set_name <- "pLI>0.9 gnomAD"
gnomadpLI <- gnomad %>%
  filter(lof.pLI > 0.9) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  pull()
  
# Add new gene set to the list
ageing_gsets <- add_gene_set(ageing_gsets, new_set_name, gnomadpLI)
length(ageing_gsets)
# Verify if the gene set has been added
print(ageing_gsets[[new_set_name]])

new_set_name <- "GTEx Up-regulated Genes"
GTExUp <- unique(GTEx$`GTExUp-regulated Genes`)

# Add new gene set to the list
ageing_gsets <- add_gene_set(ageing_gsets, new_set_name, GTExUp)
length(ageing_gsets)
# Verify if the gene set has been added
print(ageing_gsets[[new_set_name]])

new_set_name <- "GTExDown-regulated Genes"
GTExDown <- unique(GTEx$`GTExDown-regulated Genes`)

# Add new gene set to the list
ageing_gsets <- add_gene_set(ageing_gsets, new_set_name, GTExDown)
length(ageing_gsets)
# Verify if the gene set has been added
print(ageing_gsets[[new_set_name]])

# Save the new list of gene sets back to a GMT file if needed
# Function to manually write a GeneSetCollection to a GMT file
write_gmt_manual <- function(gene_sets, file_path) {
  # Open a connection to the file
  con <- file(file_path, "w")
  
  # Iterate over each gene set in the collection
  for (gene_set_name in names(gene_sets)) {
    gene_set <- gene_sets[[gene_set_name]]
    # Create a line for the GMT file
    line <- paste(gene_set_name, "", paste(geneIds(gene_set), collapse = "\t"), sep = "\t")
    # Write the line to the file
    writeLines(line, con)
  }
  
  # Close the connection to the file
  close(con)
}

# Call the function to write the GMT file
write_gmt_manual(ageing_gsets, file.path(dataDir, "agein.Hs.symbols.gmt"))

