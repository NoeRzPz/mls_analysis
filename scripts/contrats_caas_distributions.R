# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(gt)
library(ggplot2)
library(RColorBrewer)
library(UpSetR)
library(ggvenn)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load data 
# discovery data frame 
cebi_atel.caas_bygene <- read_delim(file.path(resultsDir,"/functional/cebi_atel.caas_bygene.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu.caas_bygene <- read_delim( file.path(resultsDir,"/functional/cebi_lemu.caas_bygene.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel.caas_bygene <- read_delim(file.path(resultsDir,"/functional/lemu_atel.caas_bygene.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel.caas_bygene <- read_delim(file.path(resultsDir,"/functional/cerco_atel.caas_bygene.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi.caas_bygene <- read_delim(file.path(resultsDir,"/functional/cerco_cebi.caas_bygene.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu.caas_bygene <- read_delim(file.path(resultsDir,"/functional/cerco_lemu.caas_bygene.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)



# Backgroun genes
background_genes <- read.table(file.path(resultsDir,"/caas/caastools_LQ_out/background_genes.txt"), header = FALSE, stringsAsFactors = FALSE)$V1

# spp families
spp_fam <- read_delim(file.path(dataDir,"/caas_tool_inputs/my_sp2fam_210727.tab"), 
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# CAAS distributions by Number of species found in the FG group (it excludes those ones having an indel)
contrasts <- list(cebi_atel = cebi_atel.caas_bygene, cebi_lemu = cebi_lemu.caas_bygene, lemu_atel = lemu_atel.caas_bygene, cerco_atel = cerco_atel.caas_bygene, cerco_cebi = cerco_cebi.caas_bygene, cerco_lemu = cerco_lemu.caas_bygene)

for (name in names(contrasts)) {
  # Create a one-dimensional table
  caasbygene_distribution <- table(contrasts[[name]]$CAAS_N)
  
  # Convert the table to a data frame for plotting
  df_caasbygene <- as.data.frame(caasbygene_distribution)

  # Calculate the percentages for each category
  df_caasbygene <- df_caasbygene %>%
  mutate(Percentage = Freq / sum(Freq) * 100)

  # Plot with ggplot2 using the percentages
  barp <- ggplot(data = df_caasbygene, aes(x = Var1, y = Percentage)) + 
    geom_bar(stat = "identity") +
    labs(title = paste(name, "CAAS by Gene distribution"), x = "Number CAAS/Gene", y = "Percentage") +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 100, by = 10))  # Setting y axis breaks in increments of 10
 
   # save plot
  ggsave(plot = barp, filename = file.path(resultsDir,paste0("functional/", name,".discovery/caas_barplot.png")), dpi = 300)
}

for (name in names(contrasts)) {
  # Retrieve genes with 2 or more caas
  rep_genes <- contrasts[[name]] %>%
    filter(CAAS_N >= 2) %>%
    group_by(CAAS_N) %>%
    summarise(Gene_N = n(),
              Genes = paste(Gene, collapse = ","))

  # Export table
  write.table(rep_genes, file.path(resultsDir,paste0("functional/", name,".discovery/rep_genes.tsv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}


# We select gene symbols for the ORA
genes_cebi_atel <- unique(cebi_atel.caas_bygene$Gene)
genes_cebi_lemu <- unique(cebi_lemu.caas_bygene$Gene)
genes_lemu_atel <- unique(lemu_atel.caas_bygene$Gene)
genes_cerco_atel <- unique(cerco_atel.caas_bygene$Gene)
genes_cerco_cebi <- unique(cerco_cebi.caas_bygene$Gene)
genes_cerco_lemu <- unique(cerco_lemu.caas_bygene$Gene)

# Create a list of your lists
lists <- list(CebidaeAtelidae = genes_cebi_atel, CebidaeLemuridae = genes_cebi_lemu, LemuridaeAtelidae = genes_lemu_atel, CercopithecidaeAtelidae = genes_cerco_atel, CercopithecidaeCebidae = genes_cerco_cebi, CercopithecidaeLemuridae = genes_cerco_lemu)

# Convert to binary membership data
upset_data <- fromList(lists)

# Create a named vector of colors for the sets with RColorBrewer 
set_colors <- brewer.pal(6, "Set3")
names(set_colors) <- names(lists)

# Generate the UpSet plot with customized main bar colors
upset(upset_data, sets = names(lists), sets.bar.color = set_colors, nintersects = 63)

# I check that there are genes present in all 6 gene sets, then I don't understand the way upset plot worked
#for (name in names(lists)) {
#  print(lists[[name]][which(lists[[name]] == "ABCA7")])
#}

# Use Reduce to find the intersection of all sets
shared_elements <- Reduce(intersect, lists)

# Print out the shared elements
print(shared_elements)

# Export shared genes by all gene sets
write.table(shared_elements, file.path(resultsDir,"functional/shared_genesbyall.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)


# Initialize an empty dataframe to store results
info_shared_genes <- data.frame()

# We extract CAAS in those genes for each data set, to see if they lay in the same regions and if they acumulate similar CAAS numbers
for (gene in shared_elements) {
  for (name in names(contrasts)) {
    gene_data <- contrasts[[name]] %>%
          filter(Gene == gene) %>%
      # Add a column for the dataset name
      mutate(Dataset = name) %>%
      select(Gene, Dataset, everything())

    # Append the data to the results dataframe
    info_shared_genes <- rbind(info_shared_genes, gene_data)
  }
}
# Export shared genes by all gene sets
write.table(info_shared_genes, file.path(resultsDir,"functional/info_shared_genes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Shared genes with shared CAAS are
important_genes <- c("ABCA7","ADGRG4","BAHCC1","CARD6","CCDC168","CENPJ","CMYA5","CRYBG3","CUBN","DHX37","EPB42","EPS8L3","FANCA","FAT1","FREM2", "FSIP2", "GFY","IL10RB","KNDC1","KRBA1","MXRA5","PKD1L3","PKHD1","PRRT4","PTX4","RBFA","SYNE2","SYTL2","TAS1R2","TEX15","TMC2","VWDE","WDR87","ZGRF1","ZNF469")

# Load results from ORA analysis to gather where gene of interest lay
cebi_atel.GOBP <- read_delim(file.path(resultsDir,"/functional/GOBP_cebi_atel.txt"), 
                                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu.GOBP <- read_delim( file.path(resultsDir,"/functional/GOBP_cebi_lemu.txt"), 
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel.GOBP <- read_delim(file.path(resultsDir,"/functional/GOBP_lemu_atel.txt"), 
                                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel.GOBP <- read_delim(file.path(resultsDir,"/functional/GOBP_cerco_atel.txt"), 
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi.GOBP <- read_delim(file.path(resultsDir,"/functional/GOBP_cerco_cebi.txt"), 
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu.GOBP <- read_delim(file.path(resultsDir,"/functional/GOBP_cerco_lemu.txt"), 
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

gobps <- list(CebidaeAtelidae = cebi_atel.GOBP,CebidaeLemuridae = cebi_lemu.GOBP,LemuridaeAtelidae = lemu_atel.GOBP,CercopithecidaeAtelidae = cerco_atel.GOBP,CercopithecidaeCebidae = cerco_cebi.GOBP,CercopithecidaeLemuridae = cerco_lemu.GOBP)

# Initialize the data frame for storing genes in main shared BPs
genes_mainBP <- data.frame()

for (gene in shared_elements) {
  # Initialize the data frame for the current gene
  gene_info_gobps <- data.frame()
  
  for (name in names(gobps)) {
    gene_data <- gobps[[name]] %>%
      # splits the geneID strings at each slash to create a list of character vectors, sapply() applies a function over this list to check if gene is in each vector
      filter(sapply(strsplit(as.character(geneID), "/"), function(x) gene %in% x)) %>%
      # Add a column for the dataset name
      mutate(Dataset = name) %>%
      select(Dataset, everything())
    
    # Append the data to the gene-specific dataframe
    gene_info_gobps <- rbind(gene_info_gobps, gene_data)
  }
  # Export biological processess each important gene is involved in
  write.table(gene_info_gobps, file.path(resultsDir,paste0("functional/",gene,"_info_gobps.tsv")), 
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # We filter BP statistically significant
  gene_info_gobps <- gene_info_gobps %>%
    # check if the gene is involve in the shared BPs
    filter(ID %in% c("GO:0044782","GO:0001539","GO:0060285","GO:0001578","GO:0003341","GO:1900016","GO:0048012","GO:1905246")) %>%
    # Add a column for the gene name
    mutate(gene = gene) %>%
    select(gene, everything())
  # Append the data to the gene-specific dataframe
  genes_mainBP <- rbind(genes_mainBP, gene_info_gobps)
}
# Export biological processess each important gene is involved in
write.table(genes_mainBP, file.path(resultsDir,"functional/genes_mainBP.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


# Load results from ORA GO MF analysis to gather where gene of interest lay
cebi_atel.GOMF <- read_delim(file.path(resultsDir,"/functional/GOMF_cebi_atel.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu.GOMF <- read_delim( file.path(resultsDir,"/functional/GOMF_cebi_lemu.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel.GOMF <- read_delim(file.path(resultsDir,"/functional/GOMF_lemu_atel.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel.GOMF <- read_delim(file.path(resultsDir,"/functional/GOMF_cerco_atel.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi.GOMF <- read_delim(file.path(resultsDir,"/functional/GOMF_cerco_cebi.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu.GOMF <- read_delim(file.path(resultsDir,"/functional/GOMF_cerco_lemu.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

gomfs <- list(CebidaeAtelidae = cebi_atel.GOMF,CebidaeLemuridae = cebi_lemu.GOMF,LemuridaeAtelidae = lemu_atel.GOMF,CercopithecidaeAtelidae = cerco_atel.GOMF,CercopithecidaeCebidae = cerco_cebi.GOMF,CercopithecidaeLemuridae = cerco_lemu.GOMF)

# Initialize the data frame for storing genes in main shared BPs
genes_mainMF <- data.frame()

for (gene in shared_elements) {
  # Initialize the data frame for the current gene
  gene_info_gomfs <- data.frame()
  
  for (name in names(gomfs)) {
    gene_data <- gomfs[[name]] %>%
      # splits the geneID strings at each slash to create a list of character vectors, sapply() applies a function over this list to check if gene is in each vector
      filter(sapply(strsplit(as.character(geneID), "/"), function(x) gene %in% x)) %>%
      # Add a column for the dataset name
      mutate(Dataset = name) %>%
      select(Dataset, everything())
    
    # Append the data to the gene-specific dataframe
    gene_info_gomfs <- rbind(gene_info_gomfs, gene_data)
  }
  # Export biological processess each important gene is involved in
  write.table(gene_info_gomfs, file.path(resultsDir,paste0("functional/",gene,"_info_gomfs.tsv")), 
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # We filter BP statistically significant
  gene_info_gomfs <- gene_info_gomfs %>%
    # check if the gene is involve in the shared BPs
    filter(ID %in% c("GO:0005201","GO:0019966","GO:0089720","GO:0140657","GO:0003774","GO:0005540","GO:0046978","GO:0046979","GO:1901682")) %>%
    # Add a column for the gene name
    mutate(gene = gene) %>%
    select(gene, everything())
  # Append the data to the gene-specific dataframe
  genes_mainMF <- rbind(genes_mainMF, gene_info_gomfs)
}
# Export biological processess each important gene is involved in
write.table(genes_mainMF, file.path(resultsDir,"functional/genes_mainMF.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Load results from ORA KEGG analysis to gather where gene of interest lay
cebi_atel.kegg <- read_delim(file.path(resultsDir,"/functional/kegg_cebi_atel.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu.kegg <- read_delim( file.path(resultsDir,"/functional/kegg_cebi_lemu.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel.kegg <- read_delim(file.path(resultsDir,"/functional/kegg_lemu_atel.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel.kegg <- read_delim(file.path(resultsDir,"/functional/kegg_cerco_atel.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi.kegg <- read_delim(file.path(resultsDir,"/functional/kegg_cerco_cebi.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu.kegg <- read_delim(file.path(resultsDir,"/functional/kegg_cerco_lemu.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

keggs <- list(CebidaeAtelidae = cebi_atel.kegg,CebidaeLemuridae = cebi_lemu.kegg,LemuridaeAtelidae = lemu_atel.kegg,CercopithecidaeAtelidae = cerco_atel.kegg,CercopithecidaeCebidae = cerco_cebi.kegg,CercopithecidaeLemuridae = cerco_lemu.kegg)

# Initialize the data frame for storing genes in main shared BPs
genes_mainkegg <- data.frame()

for (gene in shared_elements) {
  # Initialize the data frame for the current gene
  gene_info_keggs <- data.frame()
  
  for (name in names(keggs)) {
    gene_data <- keggs[[name]] %>%
      # splits the geneID strings at each slash to create a list of character vectors, sapply() applies a function over this list to check if gene is in each vector
      filter(sapply(strsplit(as.character(geneID), "/"), function(x) gene %in% x)) %>%
      # Add a column for the dataset name
      mutate(Dataset = name) %>%
      select(Dataset, everything())
    
    # Append the data to the gene-specific dataframe
    gene_info_keggs <- rbind(gene_info_keggs, gene_data)
  }
  # Export biological processess each important gene is involved in
  write.table(gene_info_keggs, file.path(resultsDir,paste0("functional/",gene,"_info_keggs.tsv")), 
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # We filter BP statistically significant
  gene_info_keggs <- gene_info_keggs %>%
    # check if the gene is involve in the shared BPs
    filter(ID %in% c("hsa03460","hsa04512","hsa04610")) %>%
    # Add a column for the gene name
    mutate(gene = gene) %>%
    select(gene, everything())
  # Append the data to the gene-specific dataframe
  genes_mainkegg <- rbind(genes_mainkegg, gene_info_keggs)
}
# Export biological processess each important gene is involved in
write.table(genes_mainkegg, file.path(resultsDir,"functional/genes_mainkegg.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Load results from ORA DOSE analysis to gather where gene of interest lay
cebi_atel.DO <- read_delim(file.path(resultsDir,"/functional/DO_cebi_atel.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu.DO <- read_delim( file.path(resultsDir,"/functional/DO_cebi_lemu.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel.DO <- read_delim(file.path(resultsDir,"/functional/DO_lemu_atel.txt"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel.DO <- read_delim(file.path(resultsDir,"/functional/DO_cerco_atel.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi.DO <- read_delim(file.path(resultsDir,"/functional/DO_cerco_cebi.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu.DO <- read_delim(file.path(resultsDir,"/functional/DO_cerco_lemu.txt"), 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

DOs <- list(CebidaeAtelidae = cebi_atel.DO,CebidaeLemuridae = cebi_lemu.DO,LemuridaeAtelidae = lemu_atel.DO,CercopithecidaeAtelidae = cerco_atel.DO,CercopithecidaeCebidae = cerco_cebi.DO,CercopithecidaeLemuridae = cerco_lemu.DO)

# Initialize the data frame for storing genes in main shared BPs
genes_mainDO <- data.frame()

for (gene in shared_elements) {
  # Initialize the data frame for the current gene
  gene_info_DOs <- data.frame()
  
  for (name in names(DOs)) {
    gene_data <- DOs[[name]] %>%
      # splits the geneID strings at each slash to create a list of character vectors, sapply() applies a function over this list to check if gene is in each vector
      filter(sapply(strsplit(as.character(geneID), "/"), function(x) gene %in% x)) %>%
      # Add a column for the dataset name
      mutate(Dataset = name) %>%
      select(Dataset, everything())
    
    # Append the data to the gene-specific dataframe
    gene_info_DOs <- rbind(gene_info_DOs, gene_data)
  }
  # Export biological processess each important gene is involved in
  write.table(gene_info_DOs, file.path(resultsDir,paste0("functional/",gene,"_info_DOs.tsv")), 
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # We filter DO statistically significant
  gene_info_DOs <- gene_info_DOs %>%
    # check if the gene is involve in the shared BPs
    filter(ID %in% c("DOID:0060340","DOID:0111910","DOID:12336","DOID:10140","DOID:2921","DOID:9562")) %>%
    # Add a column for the gene name
    mutate(gene = gene) %>%
    select(gene, everything())
  # Append the data to the gene-specific dataframe
  genes_mainDO <- rbind(genes_mainDO, gene_info_DOs)
}
# Export biological processess each important gene is involved in
write.table(genes_mainDO, file.path(resultsDir,"functional/genes_mainDO.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Load results from ORA DisGeNET analysis to gather where gene of interest lay
cebi_atel.DGN <- read_delim(file.path(resultsDir,"/functional/DGN_cebi_atel.txt"), 
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu.DGN <- read_delim( file.path(resultsDir,"/functional/DGN_cebi_lemu.txt"), 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel.DGN <- read_delim(file.path(resultsDir,"/functional/DGN_lemu_atel.txt"), 
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel.DGN <- read_delim(file.path(resultsDir,"/functional/DGN_cerco_atel.txt"), 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi.DGN <- read_delim(file.path(resultsDir,"/functional/DGN_cerco_cebi.txt"), 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu.DGN <- read_delim(file.path(resultsDir,"/functional/DGN_cerco_lemu.txt"), 
                            delim = "\t", escape_double = FALSE, trim_ws = TRUE)

DGNs <- list(CebidaeAtelidae = cebi_atel.DGN,CebidaeLemuridae = cebi_lemu.DGN,CercopithecidaeAtelidae = cerco_atel.DGN,CercopithecidaeCebidae = cerco_cebi.DGN,CercopithecidaeLemuridae = cerco_lemu.DGN)

# Initialize the data frame for storing genes in main shared BPs
genes_mainDGN <- data.frame()

for (gene in shared_elements) {
  # Initialize the data frame for the current gene
  gene_info_DGNs <- data.frame()
  
  for (name in names(DGNs)) {
    gene_data <- DGNs[[name]] %>%
      # splits the geneID strings at each slash to create a list of character vectors, sapply() applies a function over this list to check if gene is in each vector
      filter(sapply(strsplit(as.character(geneID), "/"), function(x) gene %in% x)) %>%
      # Add a column for the dataset name
      mutate(Dataset = name) %>%
      select(Dataset, everything())
    
    # Append the data to the gene-specific dataframe
    gene_info_DGNs <- rbind(gene_info_DGNs, gene_data)
  }
  # Export biological processess each important gene is involved in
  write.table(gene_info_DGNs, file.path(resultsDir,paste0("functional/",gene,"_info_DGNs.tsv")), 
              sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # We filter DGN statistically significant
  gene_info_DGNs <- gene_info_DGNs %>%
    # check if the gene is involve in the shared DGNs
    filter(ID %in% c("C4551493","C0021359","C0231835","C4551906")) %>%
    # Add a column for the gene name
    mutate(gene = gene) %>%
    select(gene, everything())
  # Append the data to the gene-specific dataframe
  genes_mainDGN <- rbind(genes_mainDGN, gene_info_DGNs)
}
# Export biological processess each important gene is involved in
write.table(genes_mainDGN, file.path(resultsDir,"functional/genes_mainDGN.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


# Since we only have genes shared by combinations of 3 sets we explore each option individually by Venn diagram
# Create a list of 3 lists
listCebi <- list(Cebidae.Atelidae = genes_cebi_atel, Cebidae.Lemuridae = genes_cebi_lemu, Cebidae.Cercopithecidae = genes_cerco_cebi)

ggvenn(listCebi, stroke_size = 0.5, set_name_size = 3)


# Function to get all combinations of 6 sets from the list of 6 sets
combinations_of_five <- combn(lists, 6, simplify = FALSE)
# Function to intersect each combination of five sets
intersections_of_five <- lapply(combinations_of_five, function(sets) Reduce(intersect, sets))
# Flatten the list of intersections to a single vector
all_intersections <- unlist(intersections_of_five, recursive = FALSE)
# Create a table of counts for each element
element_counts <- table(unlist(all_intersections))
element_counts


df <- data.frame(
  element = unlist(lists),
  set = rep(names(lists), times = lapply(lists, length))
)

# Count occurrences of each element
element_counts <- df %>%
  group_by(element) %>%
  summarize(count = n_distinct(set))

# Filter elements present in exactly five sets
elements_in_five_sets <- element_counts %>%
  filter(count == 5) %>%
  pull(element)

# Print the elements
unique(elements_in_five_sets)
