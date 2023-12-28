# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(gt)
library(ggplot2)
library(RColorBrewer)
library(UpSetR)
library(ComplexUpset)
library(ggvenn)
library(cowplot)
library(patchwork)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load data from  family contrasts 
cebi_atel <- read_delim(file.path(resultsDir,"/functional/cebi_atel_byFGN.txt"),col_names = TRUE ,
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu <- read_delim(file.path(resultsDir,"/functional/cebi_lemu_byFGN.txt"),col_names = TRUE ,
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel <- read_delim(file.path(resultsDir,"/functional/lemu_atel_byFGN.txt"),col_names = TRUE ,
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel <- read_delim(file.path(resultsDir,"/functional/cerco_atel_byFGN.txt"),col_names = TRUE ,
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi <- read_delim(file.path(resultsDir,"/functional/cerco_cebi_byFGN.txt"),col_names = TRUE ,
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu <- read_delim(file.path(resultsDir,"/functional/cerco_lemu_byFGN.txt"),col_names = TRUE ,
                         delim = "\t", escape_double = FALSE, trim_ws = TRUE)

fam4 <- read_delim(file.path(resultsDir,"/functional/4fam_byFGN.txt"),col_names = TRUE ,
                   delim = "\t", escape_double = FALSE, trim_ws = TRUE)


# Obtain gene lists without duplicates
genes_cebi_atel <- unique(cebi_atel$Gene)

genes_cebi_lemu <- unique(cebi_lemu$Gene)

genes_lemu_atel <- unique(lemu_atel$Gene)

genes_cerco_atel <- unique(cerco_atel$Gene)

genes_cerco_cebi <- unique(cerco_cebi$Gene)

genes_cerco_lemu <- unique(cerco_lemu$Gene)

genes_4fam <- unique(fam4$Gene)

# Create a list of your lists
lists <- list(Cebidae.Atelidae = genes_cebi_atel, Cebidae.Lemuridae = genes_cebi_lemu, Lemuridae.Atelidae = genes_lemu_atel, Cercopithecidae.Atelidae = genes_cerco_atel, Cercopithecidae.Cebidae = genes_cerco_cebi, Cercopithecidae.Lemuridae = genes_cerco_lemu, All4fam = genes_4fam)

# Convert to binary membership data
upset_data <- fromList(lists)

# Create a named vector of colors for the sets with RColorBrewer 
set_colors <- brewer.pal(7, "Set3")
names(set_colors) <- names(lists)

# Calculate gene overlap between family contrasts
# Generate the UpSet plot with customized colors
u7 <- upset(upset_data, sets = names(lists), sets.bar.color = set_colors, nintersects = 80)

# Open a PNG device
png(filename = file.path(resultsDir, "/descriptive/genes_upsetplot_7contrasts.png"), width = 13, height = 7, units = 'in', res = 500)
# Draw the plot on the open device
print(u7)
# Turn off the device to close and save the file
dev.off()

# Create the upset plot with ComplexUpset
u7 <- ComplexUpset::upset(upset_data, intersect = names(lists), width_ratio = 0.1, n_intersections = 80,name = '',
                          base_annotations = list('Intersection size' = (
                            ComplexUpset::intersection_size(
                              text = list(vjust = -0.1,color = 'black',hjust = -0.1,angle = 45, size = 3)) +
                              theme_classic() + theme( axis.title.x = element_blank(),axis.text.x = element_blank()))))


# Create a list to know total CAAS discovered with all contrasts
caas_cebi_atel <- cebi_atel %>%
  mutate(caas = paste(Gene, Position, sep = "_")) %>%
  dplyr::select(Gene,caas)

caas_cebi_lemu <- cebi_lemu %>%
  mutate(caas = paste(Gene, Position, sep = "_")) %>%
  dplyr::select(Gene,caas)

caas_lemu_atel <- lemu_atel %>%
  mutate(caas = paste(Gene, Position, sep = "_")) %>%
  dplyr::select(Gene,caas)

caas_cerco_atel <- cerco_atel %>%
  mutate(caas = paste(Gene, Position, sep = "_")) %>%
  dplyr::select(Gene,caas)

caas_cerco_cebi <- cerco_cebi %>%
  mutate(caas = paste(Gene, Position, sep = "_")) %>%
  dplyr::select(Gene,caas)

caas_cerco_lemu <- cerco_lemu %>%
  mutate(caas = paste(Gene, Position, sep = "_")) %>%
  dplyr::select(Gene,caas)

caas_fam4 <- fam4 %>%
  mutate(caas = paste(Gene, Position, sep = "_")) %>%
  dplyr::select(Gene,caas)

#joing all, and remove the duplicated to know the number of unique cass
all7contrasts <-  rbind(caas_cebi_atel, caas_cebi_lemu, caas_lemu_atel, caas_cerco_atel, caas_cerco_cebi, caas_cerco_lemu, caas_fam4)

#Unique CAAS
length(unique(sort(all7contrasts$caas)))
length(unique(sort(caas_cebi_lemu$caas)))

# unique genes
length(unique(sort(all7contrasts$Gene)))
length(unique(sort(caas_cebi_lemu$Gene)))


lists <- list(all4fam = caas_fam4$caas, Cebidae.Atelidae = caas_cebi_atel$caas, Cebidae.Lemuridae = caas_cebi_lemu$caas, Lemuridae.Atelidae = caas_lemu_atel$caas, Cercopithecidae.Atelidae = caas_cerco_atel$caas, Cercopithecidae.Cebidae = caas_cerco_cebi$caas, Cercopithecidae.Lemuridae = caas_cerco_lemu$caas)

# Convert to binary membership data
upset_data <- fromList(lists)

# Create a named vector of colors for the sets with RColorBrewer 
set_colors <- brewer.pal(7, "Set3")
names(set_colors) <- names(lists)

# Generate the UpSet plot with customized colors
c7 <- upset(upset_data, sets = names(lists), sets.bar.color = set_colors)

# Open a PNG device
png(filename = file.path(resultsDir, "/descriptive/caas_upsetplot_7contrasts.png"), width = 13, height = 7, units = 'in', res = 500)
# Draw the plot on the open device
print(c7)
# Turn off the device to close and save the file
dev.off()

# Create the upset plot with ComplexUpset
c7 <- ComplexUpset::upset(upset_data, intersect = names(lists), width_ratio = 0.2, n_intersections = 80,name = '',
                          base_annotations = list('Intersection size' = (
                            ComplexUpset::intersection_size(
                              text = list(vjust = -0.1,color = 'black',hjust = -0.1,angle = 45)) +
                              theme_classic() + theme( axis.title.x = element_blank(),axis.text.x = element_blank()))))
# Use Reduce to find the intersection of all sets
shared_elements <- Reduce(intersect, lists)

# Print out the shared elements
shared_elements

# If c7 and u7 are lists of ggplot objects, use wrap_plots to combine each list into a single plot
plot_c7 <- plot_grid(c7, ncol = 1) 
plot_u7 <- plot_grid(u7, ncol = 1) 

# Combine the two annotated plot objects
combined_plot <- plot_c7 / plot_u7  + 
  plot_annotation(tag_levels = 'A')
combined_plot
# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"/descriptive/ggplot_upsetplots_7contrasts.png"),  dpi = 600)

# Create figure 2
# Read the images
image1 <- magick::image_read(file.path(resultsDir,"/descriptive/caas_upsetplot_7contrasts.png"))
image2 <- magick::image_read(file.path(resultsDir,"/descriptive/genes_upsetplot_7contrasts.png"))

# Annotate each image with letters A, B, and C respectively
image1 <- ggdraw() + draw_image(image1)
image2 <- ggdraw() + draw_image(image2)

combined_plot <- (image1 / image2) + 
  plot_annotation(tag_levels = 'A')
combined_plot
# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"/descriptive/combined_upsetplots_7contrasts.png"),  dpi = 600)


# Use Reduce to find the intersection of all sets
shared_elements <- Reduce(intersect, lists)

# Print out the shared elements
shared_elements


df <- data.frame(
  element = unlist(lists),
  set = rep(names(lists), times = lapply(lists, length))
)

# Count occurrences of each element
element_counts <- df %>%
  group_by(element) %>%
  summarize(count = n_distinct(set))

# Filter elements present in exactly five sets
elements_in_7_sets <- element_counts %>%
  filter(count == 7) %>%
  pull(element)

# Print the elements
print(elements_in_7_sets)


# Create a list of your lists
lists <- list(Cebidae.Atelidae = genes_cebi_atel, Cebidae.Lemuridae = genes_cebi_lemu, Lemuridae.Atelidae = genes_lemu_atel, Cercopithecidae.Atelidae = genes_cerco_atel, Cercopithecidae.Cebidae = genes_cerco_cebi, Cercopithecidae.Lemuridae = genes_cerco_lemu)

# Convert to binary membership data
upset_data <- fromList(lists)

# Create a named vector of colors for the sets with RColorBrewer 
set_colors <- brewer.pal(6, "Set3")
names(set_colors) <- names(lists)

# Generate the UpSet plot with customized main bar colors
u6 <- upset(upset_data, sets = names(lists), sets.bar.color = set_colors, nintersects = 63)

# Open a PNG device
png(filename = file.path(resultsDir, "/descriptive/genes_upsetplot_6contrasts.png"), width = 13, height = 7, units = 'in', res = 500)
# Draw the plot on the open device
print(u6)
# Turn off the device to close and save the file
dev.off()

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
cebi_atel.GOBP <- read_delim(file.path(resultsDir,"/functional/ORA_results/GOBP_cebi_atel.txt"), 
                                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu.GOBP <- read_delim( file.path(resultsDir,"/functional/ORA_results/GOBP_cebi_lemu.txt"), 
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel.GOBP <- read_delim(file.path(resultsDir,"/functional/ORA_results/GOBP_lemu_atel.txt"), 
                                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel.GOBP <- read_delim(file.path(resultsDir,"/functional/ORA_results/GOBP_cerco_atel.txt"), 
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi.GOBP <- read_delim(file.path(resultsDir,"/functional/ORA_results/GOBP_cerco_cebi.txt"), 
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu.GOBP <- read_delim(file.path(resultsDir,"/functional/ORA_results/GOBP_cerco_lemu.txt"), 
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
