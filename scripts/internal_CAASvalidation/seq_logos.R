# Sequences Logos for validated genes
# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(ggplot2)
library(ggplot2)
library(ggmsa)
library(Biostrings)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

mammals <- read_tsv(file.path(resultsDir,"/functional/external_validation/labels.txt"), col_names = T)
# Get the list of files that match the pattern
file_names <- list.files(path = file.path(resultsDir, "/functional/external_validation/MultipleCodonAlignments"),
                         pattern = paste0("ENST.*", "SH3PXD2B", ".fasta.gz$"), full.names = TRUE)

# Check if files are found
if (length(file_names) == 0) {
  stop("No files found for the specified pattern")
}

# Read the first file (or loop through all files if necessary)
myMsaObject <- Biostrings::readDNAMultipleAlignment(file_names[1], format = "fasta")
# Function to translate each sequence in the alignment
translate_alignment <- function(dna_alignment) {
  # Convert the DNAMultipleAlignment object to a list of DNAString objects
  dna_seq_list <- as(dna_alignment, "List")
  
  # Use lapply to translate each sequence
  lapply(dna_seq_list, function(seq) {
    # Perform the translation
    aa_seq <- translate(seq)
    # Return the translated sequence, removing any stop codons
    #removeStops(aa_seq)
  })
}

# Translate the DNA alignment
protein_alignment_list <- translate_alignment(myMsaObject)

# Optional: Convert the list back to a MultipleSequenceAlignment object
protein_alignment <- MultipleAlignment(protein_alignment_list)

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[1:10]

# Generate plot
ggmsa(myMsaObject, start = 730, end = 740, char_width = 0.5, seq_name = TRUE) 



# We define all extreme LQ species detected by family quantile criteria
# Ordered from top to bottom LQ
spp <- c('Sapajus_apella','Ateles_geoffroyi','Macaca_fascicularis','Macaca_fuscata','Macaca_mulatta','Macaca_nemestrina', 'Eulemur_mongoz', 'Homo_sapiens',  'Cheirogaleus_medius','Saguinus_imperator','Alouatta_palliata','Nasalis_larvatus','Trachypithecus_francoisi','Semnopithecus_entellus','Pygathrix_nemaeus','Prolemur_simus') 

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "SH3PXD2B", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 730, end = 740, char_width = 0.5, seq_name = TRUE) 


# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "PRR14L", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 1437, end = 1448, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "JPH3", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 580, end = 591, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "SEL1L2", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 10, end = 25, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()


# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "HYDIN", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Subset and reorder the MSA object (reorder w/o eulemur mongoz and Prolemur simus that seem to be ausent)
ordered <-  c('Sapajus_apella','Macaca_fascicularis','Macaca_fuscata','Macaca_mulatta','Macaca_nemestrina', 'Eulemur_mongoz', 'Homo_sapiens',  'Cheirogaleus_medius','Saguinus_imperator','Nasalis_larvatus','Trachypithecus_francoisi','Semnopithecus_entellus','Pygathrix_nemaeus','Prolemur_simus') 
myMsaObject@unmasked <- myMsaObject@unmasked[ordered]

# Generate plot
ggmsa(myMsaObject, start = 2525, end = 2536, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "NCAPH2", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 287, end = 298, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "PLB1", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 1344, end = 1355, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()


# Load genes shared in all 2 family contrasts that passed internal validation 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "DNAH8", ".Homo_sapiens.filter2.phy")), format = "phylip")
# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]

# Subset and reorder the MSA object (reorder w/o eulemur mongoz and Prolemur simus that seem to be ausent)
ordered <- c('Sapajus_apella','Ateles_geoffroyi','Macaca_fascicularis','Macaca_fuscata','Macaca_mulatta','Macaca_nemestrina', 'Homo_sapiens',  'Cheirogaleus_medius','Saguinus_imperator','Alouatta_palliata','Nasalis_larvatus','Trachypithecus_francoisi','Semnopithecus_entellus','Pygathrix_nemaeus') 
myMsaObject@unmasked <- myMsaObject@unmasked[ordered]
# Generate plot
ggmsa(myMsaObject, start = 2863, end = 2873, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "SRRM5", ".Homo_sapiens.filter2.phy")), format = "phylip")
# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]

# Subset and reorder the MSA object (reorder w/o eulemur mongoz and Prolemur simus that seem to be ausent)
ordered <- c('Sapajus_apella','Ateles_geoffroyi','Macaca_fascicularis','Macaca_fuscata','Macaca_mulatta','Macaca_nemestrina', 'Homo_sapiens',  'Cheirogaleus_medius','Saguinus_imperator','Alouatta_palliata','Nasalis_larvatus','Trachypithecus_francoisi','Semnopithecus_entellus','Pygathrix_nemaeus') 
myMsaObject@unmasked <- myMsaObject@unmasked[ordered]
# Generate plot
ggmsa(myMsaObject, start = 430, end = 441, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()


myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "CUZD1", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 340, end = 351, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "METTL22", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 217, end = 227, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "NSRP1", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 464, end = 475, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "ASB1", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 152, end = 163, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "GBP1", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 324, end = 335, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "CACNA1H", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Reorder the MSA object 
myMsaObject@unmasked <- myMsaObject@unmasked[spp]

# Generate plot
ggmsa(myMsaObject, start = 1849, end = 1860, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "DDX54", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
# Subset and reorder the MSA object (reorder w/o Sapajus apella)
ordered <- c('Ateles_geoffroyi','Macaca_fascicularis','Macaca_fuscata','Macaca_mulatta','Macaca_nemestrina', 'Homo_sapiens',  'Cheirogaleus_medius','Saguinus_imperator','Alouatta_palliata','Nasalis_larvatus','Trachypithecus_francoisi','Semnopithecus_entellus','Pygathrix_nemaeus') 
myMsaObject@unmasked <- myMsaObject@unmasked[ordered]

# Generate plot
ggmsa(myMsaObject, start = 45, end = 56, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# Load data 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "ALPK1", ".Homo_sapiens.filter2.phy")), format = "phylip")
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "GNRH2", ".Homo_sapiens.filter2.phy")), format = "phylip")

# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]

ggmsa(myMsaObject, start = 5, end = 15, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "GNRH2", ".Homo_sapiens.filter2.phy")), format = "phylip")
# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
ggmsa(myMsaObject, start = 1, end = 15, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "GNRH2", ".Homo_sapiens.filter2.phy")), format = "phylip")
# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
ggmsa(myMsaObject, start = 1, end = 15, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()
 
myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "HSD3B1", ".Homo_sapiens.filter2.phy")), format = "phylip")
# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
ggmsa(myMsaObject, start = 10, end = 23, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "IQCJ-SCHIP1", ".Homo_sapiens.filter2.phy")), format = "phylip")
# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
ggmsa(myMsaObject, start = 85, end = 100, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "MAGEA11", ".Homo_sapiens.filter2.phy")), format = "phylip")
# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
ggmsa(myMsaObject, start = 200, end = 215, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

myMsaObject <- Biostrings::readAAMultipleAlignment(file.path(dataDir, paste0("/palignments.231026/", "ZNF107", ".Pygathrix_nemaeus.filter2.phy")), format = "phylip")
# Subset the MSA object
myMsaObject@unmasked <- myMsaObject@unmasked[names(myMsaObject@unmasked) %in% spp]
ggmsa(myMsaObject, start = 220, end = 231, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()
ggmsa(myMsaObject, start = 310, end = 321, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()

# To save the plot
ggsave("alignment_plot.pdf", plot = p, width = 11, height = 8.5, dpi = 300)
