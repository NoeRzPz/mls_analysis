# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(gt)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler) 
library(DOSE)
library(msigdbr)
library(pheatmap)
#library(ape)
library(ReactomePA)
library(cowplot)
library(patchwork)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load data 
# discovery data frame 
ceat_leche.discovery <- read_delim(file.path(resultsDir,"/caas/LQ_cebi_ateli_vs_lemu_cheiro/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_atel.discovery <- read_delim(file.path(resultsDir,"/caas/LQ_cebidae_atelidae/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cebi_lemu.discovery <- read_delim(file.path(resultsDir,"/caas/LQ_cebidae_lemuridae/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

lemu_atel.discovery <- read_delim(file.path(resultsDir,"/caas/LQ_lemuridae_atelidae/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_atel.discovery <- read_delim(file.path(resultsDir,"/caas/LQ_cercopi_atelidae/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_cebi.discovery <- read_delim(file.path(resultsDir,"/caas/LQ_cercopi_cebidae/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

cerco_lemu.discovery <- read_delim(file.path(resultsDir,"/caas/LQ_cercopi_lemuridae/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Read in the gene lengths data
gene_lengths <- read_delim(file.path(resultsDir,"/functional/gene_lengths.tsv"),  delim = "\t")

# Backgroun genes
background_genes <- read.table(file.path(resultsDir,"/caas/caastools_LQ_out/background_genes.txt"), header = FALSE, stringsAsFactors = FALSE)$V1

# spp families
spp_fam <- read_delim(file.path(dataDir,"/caas_tool_inputs/my_sp2fam_210727.tab"), 
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# 70 genes present in all 2 family contrasts
genes_shared_all2fam <- read_tsv(file.path(resultsDir,"functional/shared_genesbyall.tsv"), col_names = F)$X1

# CAAS distributions by Number of species found in the FG group (it excludes those ones having an indel)
contrasts <- list(ceat_leche.discovery = ceat_leche.discovery, cebi_atel.discovery = cebi_atel.discovery, cebi_lemu.discovery = cebi_lemu.discovery, lemu_atel.discovery = lemu_atel.discovery, cerco_atel.discovery = cerco_atel.discovery, cerco_cebi.discovery = cerco_cebi.discovery, cerco_lemu.discovery = cerco_lemu.discovery)

for (name in names(contrasts)) {
  
  freq_table <- contrasts[[name]] %>%
  group_by(FFGN) %>%
  summarise(CAAS_N = n()) %>%
  arrange(FFGN)
  
  # Define the path for the new directory
  dir_path <- file.path(resultsDir,paste0("/functional/", name,"/"))
  
  # Check if the directory already exists
  if (!dir.exists(dir_path)) {
    # Create the directory since it does not exist
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Export table
  write.table(freq_table , file.path(resultsDir,paste0("/functional/", name,"/CAAS_byFFGN.txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}

for (name in names(contrasts)) {
  freq_table <-  contrasts[[name]] %>%
  filter(Pvalue <= 0.05) %>%
  count(FBGN, FFGN) %>%
  spread(key = FFGN, value = n, fill = 0)
  
  # Formatting the table using the gt package
  formatted_table <- freq_table %>%
    gt() %>%
    tab_header(title = "FFGN") %>%
    cols_label( FBGN = "FBGN")
  
  # Export table
  write.table(formatted_table, file.path(resultsDir,paste0("/functional/", name,"/CAAS_byFFGN005.txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}
  

for (name in names(contrasts)) {
  cat(paste0(name,"\n"))
  print(dim(contrasts[[name]])[1])
}

for (name in names(contrasts)) {
  # Explore what are the species contributing more to each case
  freq_spp <- contrasts[[name]] %>%
    group_by(FFGN) %>%
    # Separate the species into individual rows
    separate_rows(FFG, sep = ",") %>%
    mutate(FFG = str_replace(FFG, "_", " ")) %>%
    # Count the occurrences of each species
    count(FFG)  %>% 
    ungroup() %>% # Remove group_by, as it's not needed anymore
    # Pivot to wide format
    pivot_wider(
      names_from = FFG, 
      values_from = n,
      values_fill = list(n = 0)  # Fill missing counts with 0
    ) %>%
    # Arrange by FFGN 
    arrange(FFGN)
  
  # Convert the dataframe to a matrix, excluding the FFGN column
  heatmap_data <- as.matrix(freq_spp[,-1])
  
  # Add the FFGN as row names to the matrix
  rownames(heatmap_data) <- freq_spp$FFGN

  # Ensure the species names match the column names of your heatmap data
  heatmap_data_colnames <- colnames(heatmap_data)

  # Create the annotation dataframe with family info
  annotation <- spp_fam %>%
    mutate(Family = str_replace(Family, "\\w+_", ""),
           Species = str_replace(Species, "_", " ")) %>%
    filter(Species %in% heatmap_data_colnames) %>%
    column_to_rownames(var = "Species")

  # Assuming 'unique_families' is a vector of all unique family names from your 'species_families' data
  unique_families <- sort(unique(annotation$Family))
  
  # I define my own palette
  famcol <- c("#f8766d",  "#a3a500", "#08A045","#9590ff", "#d89000")
  
  # Create a color vector where each family is assigned an hexadecimal value
  family_colors <- setNames(famcol[seq_along(unique_families)], unique_families)

  # Create a list for annotation_colors
  annotation_colors <- list(Family = family_colors)
  
  # Create heatmap. Escale by number of foreground species found where CAAS were found
  p <- pheatmap(heatmap_data, scale = "row", cluster_rows = FALSE,  cluster_cols = TRUE,
         annotation_col = annotation,
         annotation_colors = annotation_colors,
         angle_col = 315, # To show column names with an angle
         fontsize_row = 12) # To show rows with bigger font size

  # Use traditional graphics device to save the heatmap
  png(file.path(resultsDir,paste0("functional/", name,"/FGheatmap.png")), width = 6, height = 4, units = 'in', res = 300)
  print(p)
  dev.off()

}
for (name in names(contrasts)) {
  # Explore what are the species contributing more to each case
  freq_spp <- contrasts[[name]] %>%
    group_by(FBGN) %>%
    # Separate the species into individual rows
    separate_rows(FBG, sep = ",") %>%
    mutate(FBG = str_replace(FBG, "_", " ")) %>%
    # Count the occurrences of each species
    count(FBG)  %>% 
    ungroup() %>% # Remove group_by, as it's not needed anymore
    # Pivot to wide format
    pivot_wider(
      names_from = FBG, 
      values_from = n,
      values_fill = list(n = 0)  # Fill missing counts with 0
    ) %>%
    # Arrange by FBGN 
    arrange(FBGN)
  
  # Convert the dataframe to a matrix, excluding the FBGN column
  heatmap_data <- as.matrix(freq_spp[,-1])
  
  # Add the FBGN as row names to the matrix
  rownames(heatmap_data) <- freq_spp$FBGN
  
  # Ensure the species names match the column names of your heatmap data
  heatmap_data_colnames <- colnames(heatmap_data)
  
  # Create the annotation dataframe with family info
  annotation <- spp_fam %>%
    mutate(Family = str_replace(Family, "\\w+_", ""),
           Species = str_replace(Species, "_", " ")) %>%
    filter(Species %in% heatmap_data_colnames) %>%
    column_to_rownames(var = "Species")
  
  # Assuming 'unique_families' is a vector of all unique family names from your 'species_families' data
  unique_families <- sort(unique(annotation$Family))
  
  # I define my own palette
  famcol <- c("#f8766d",  "#a3a500", "#08A045","#9590ff", "#d89000")
  
  # Create a color vector where each family is assigned an hexadecimal value
  family_colors <- setNames(famcol[seq_along(unique_families)], unique_families)
  
  # Create a list for annotation_colors
  annotation_colors <- list(Family = family_colors)
  
  # Create heatmap. Escale by number of foreground species found where CAAS were found
  p <- pheatmap(heatmap_data, scale = "row", cluster_rows = FALSE,  cluster_cols = TRUE,
                annotation_col = annotation,
                annotation_colors = annotation_colors,
                angle_col = 315, # To show column names with an angle
                fontsize_row = 12) # To show rows with bigger font size
  
  # Use traditional graphics device to save the heatmap
  png(file.path(resultsDir,paste0("functional/", name,"/BGheatmap.png")), width = 6, height = 4, units = 'in', res = 300)
  print(p)
  dev.off()
  
}

for (name in names(contrasts)) {
  # Plot histogram of hipergeometric test pvalues
  pva_histo <- ggplot(data = contrasts[[name]], aes(x = Pvalue)) + 
    geom_histogram(binwidth = 0.02, fill = "blue", color = "black", alpha = 0.7) + 
    theme_classic() +
    labs(title = paste("Histogram of P-values", name), x = "P-value", y = "Frequency")

  # save plot
  ggsave(plot = pva_histo, filename = file.path(resultsDir,paste0("functional/", name,"/pval_histo.png")), dpi = 300)
}

for (name in names(contrasts)) {
  # Scatter plot
  fg_bg <- contrasts[[name]] %>% 
    ggplot(aes(x = FBGN, y = FFGN, color = FFGN)) +
    geom_jitter() +
    labs( x = 'FBGN',
          y = 'FFGN') +
    theme_classic() + 
    theme(
      legend.position = "bottom",
      text = element_text(size = 16))

  # save plot
  ggsave(plot = fg_bg, filename = file.path(resultsDir,paste0("functional/", name,"/fg_bgscatter.png")), dpi = 300)
}

for (name in names(contrasts)) {
  # CAAS distributions by Number of species found in the FG group. Count occurrences
  freq_table <-  contrasts[[name]] %>%
    count(FBGN, FFGN) %>%
    spread(key = FFGN, value = n, fill = 0)

  # Formatting the table using the gt package
  formatted_table <- freq_table %>%
    gt() %>%
    tab_header(title = "FFGN") %>%
    cols_label(
      FBGN = "FBGN")

  # Export table, buscar como era con gt
  write.table(formatted_table, file.path(resultsDir,paste0("functional/", name,"/CAAS_byFFGN.txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}

for (name in names(contrasts)) {
  # Count occurrences gap
  freq_table <-  contrasts[[name]] %>%
    count(GBG, GFG) %>%
    spread(key = GFG, value = n, fill = 0)
  
  # Formatting the table using the gt package
  formatted_table <- freq_table %>%
    gt() %>%
    tab_header(title = "GFG") %>%
    cols_label(
      GBG = "GBG")
  
  # Export table, buscar como era con gt
  write.table(formatted_table, file.path(resultsDir,paste0("functional/", name,"/CAAS_byGaps.txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}

# First we performed analysis w/o filtering by nominal pvalue
# Pattern distribution in nominal significant results
freq_table <- caas.discovery %>%
  group_by(Pattern) %>%
  summarise(Frequency = n(),
            Genes = paste(unique(Gene), collapse = ",")) %>%
  arrange(Pattern)

# We filter each discovery dataframe contrast by number of FG and BG species in which a CAAS was detected

ceat_leche_byFGN <- ceat_leche.discovery %>%
  filter(FFGN >= 3  & FBGN == 3) 

cebi_atel_byFGN <- cebi_atel.discovery %>%
  filter(FFGN == 2  & FBGN == 2) 

cebi_lemu_byFGN <- cebi_lemu.discovery %>%
  filter(FFGN == 2  & FBGN == 2)
write.table(cebi_lemu_byFGN, file.path(resultsDir,"/functional/cebi_lemu_discovery.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

lemu_atel_byFGN <- lemu_atel.discovery %>%
  filter(FFGN == 2  & FBGN == 2) 

cerco_atel_byFGN <- cerco_atel.discovery %>%
  filter(FFGN == 4  & FBGN == 4) 

cerco_cebi_byFGN <- cerco_cebi.discovery %>%
  filter(FFGN == 4  & FBGN == 4) 

cerco_lemu_byFGN <- cerco_lemu.discovery %>%
  filter(FFGN == 4  & FBGN == 4) 

# Count caas and gaps present
ceat.leche_byFGN <- ceat_leche_byFGN %>%
  summarise(CAAS_N = length(Substitution), 
            GFG_N = sum(GFG),
            GBG_N = sum(GBG),
            Pval_0.05 = sum(Pvalue <= 0.05))

cebi.atel_byFGN <- cebi_atel_byFGN %>%
  summarise(CAAS_N = length(Substitution), 
            GFG_N = sum(GFG),
            GBG_N = sum(GBG),
            Pval_0.05 = sum(Pvalue <= 0.05))

cebi.lemu_byFGN <- cebi_lemu_byFGN %>%
  summarise(CAAS_N = length(Substitution), 
            GFG_N = sum(GFG),
            GBG_N = sum(GBG),
            Pval_0.05 = sum(Pvalue <= 0.05))

ceat_leche_byFGN <- ceat_leche_byFGN %>%
  filter(Pvalue <= 0.05)
# Export table
write.table(ceat_leche_byFGN, file.path(resultsDir,"/functional/ceat_leche_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cebi_atel_byFGN <- cebi_atel_byFGN %>%
  filter(Pvalue <= 0.05)
# Export table
write.table(cebi_atel_byFGN, file.path(resultsDir,"/functional/cebi_atel_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cebi_lemu_byFGN <- cebi_lemu_byFGN %>%
  filter(Pvalue <= 0.05)
# Export table
write.table(cebi_lemu_byFGN, file.path(resultsDir,"/functional/cebi_lemu_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

lemu_atel_byFGN <- lemu_atel_byFGN %>%
  filter(Pvalue <= 0.05)
# Export table
write.table(lemu_atel_byFGN, file.path(resultsDir,"/functional/lemu_atel_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cerco_atel_byFGN <- cerco_atel_byFGN %>%
  filter(Pvalue <= 0.05)
# Export table
write.table(cerco_atel_byFGN, file.path(resultsDir,"/functional/cerco_atel_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cerco_cebi_byFGN <- cerco_cebi_byFGN %>%
  filter(Pvalue <= 0.05)
# Export table
write.table(cerco_cebi_byFGN, file.path(resultsDir,"/functional/cerco_cebi_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cerco_lemu_byFGN <- cerco_lemu_byFGN %>%
  filter(Pvalue <= 0.05)
# Export table
write.table(cerco_lemu_byFGN, file.path(resultsDir,"/functional/cerco_lemu_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# We identify where are the gaps
ceat_leche_gaps <- ceat_leche_byFGN %>%
  filter(GFG != 0 | GBG != 0) 
ceat_leche_gaps

cebi_atel_gaps <- cebi_atel_byFGN %>% # There are no gaps
  filter(GFG != 0 | GBG != 0) 
cebi_atel_gaps

cebi_lemu_gaps <- cebi_lemu_byFGN %>% # There are no gaps
  filter(GFG != 0 | GBG != 0) 
cebi_lemu_gaps

lemu_atel_gaps <- lemu_atel_byFGN %>% # There are no gaps
  filter(GFG != 0 | GBG != 0) 
lemu_atel_gaps

cerco_atel_gaps <- cerco_atel_byFGN %>% # There are no gaps
  filter(GFG != 0 | GBG != 0) 
cerco_atel_gaps

cerco_cebi_gaps <- cerco_cebi_byFGN %>% # There are no gaps
  filter(GFG != 0 | GBG != 0) 
cerco_cebi_gaps

cerco_lemu_gaps <- cerco_lemu_byFGN %>% # There are no gaps
  filter(GFG != 0 | GBG != 0) 
cerco_lemu_gaps

# We explore the distribution of CAAS by gene, to detect if there are genes that accumulate more CAAS
cebi_atel.caas_bygene <- cebi_atel_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            Position = paste(Position, collapse = ";"),
            Substitution = paste(Substitution, collapse = ";"))
# Export table
write.table(cebi_atel.caas_bygene, file.path(resultsDir,"/functional/cebi_atel.caas_bygene.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cebi_lemu.caas_bygene <- cebi_lemu_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            Position = paste(Position, collapse = ";"),
            Substitution = paste(Substitution, collapse = ";"))
# Export table
write.table(cebi_lemu.caas_bygene, file.path(resultsDir,"/functional/cebi_lemu.caas_bygene.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

lemu_atel.caas_bygene <- lemu_atel_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            Position = paste(Position, collapse = ";"),
            Substitution = paste(Substitution, collapse = ";"))
# Export table
write.table(lemu_atel.caas_bygene, file.path(resultsDir,"/functional/lemu_atel.caas_bygene.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cerco_atel.caas_bygene <- cerco_atel_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            Position = paste(Position, collapse = ";"),
            Substitution = paste(Substitution, collapse = ";"))
# Export table
write.table(cerco_atel.caas_bygene, file.path(resultsDir,"/functional/cerco_atel.caas_bygene.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cerco_cebi.caas_bygene <- cerco_cebi_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            Position = paste(Position, collapse = ";"),
            Substitution = paste(Substitution, collapse = ";"))
# Export table
write.table(cerco_cebi.caas_bygene, file.path(resultsDir,"/functional/cerco_cebi.caas_bygene.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

cerco_lemu.caas_bygene <- cerco_lemu_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            Position = paste(Position, collapse = ";"),
            Substitution = paste(Substitution, collapse = ";"))
# Export table
write.table(cerco_lemu.caas_bygene, file.path(resultsDir,"/functional/cerco_lemu.caas_bygene.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


ceat_leche_pat <- ceat_leche_byFGN %>%
  group_by(Pattern) %>%
  summarise(CAAS = n()) %>%
  arrange(Pattern)
ceat_leche_pat

cebi_atel_pat <- cebi_atel_byFGN %>%
  group_by(Pattern) %>%
  summarise(CAAS = n()) %>%
  arrange(Pattern)
cebi_atel_pat

cebi_lemu_pat <- cebi_lemu_byFGN %>%
  group_by(Pattern) %>%
  summarise(CAAS = n()) %>%
  arrange(Pattern)
cebi_lemu_pat

lemu_atel_pat <- lemu_atel_byFGN %>%
  group_by(Pattern) %>%
  summarise(CAAS = n()) %>%
  arrange(Pattern)
lemu_atel_pat

cerco_atel_pat <- cerco_atel_byFGN %>%
  group_by(Pattern) %>%
  summarise(CAAS = n()) %>%
  arrange(Pattern)
cerco_atel_pat

cerco_cebi_pat <- cerco_cebi_byFGN %>%
  group_by(Pattern) %>%
  summarise(CAAS = n()) %>%
  arrange(Pattern)
cerco_cebi_pat

cerco_lemu_pat <- cerco_lemu_byFGN %>%
  group_by(Pattern) %>%
  summarise(CAAS = n()) %>%
  arrange(Pattern)
cerco_lemu_pat

write.table(df_byFGN , file.path(resultsDir,"/functional/df_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# We select gene symbols for the ORA
genes_ceat_leche <- unique(ceat_leche_byFGN$Gene)
# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(genes_ceat_leche, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

genes_cebi_atel <- unique(cebi_atel_byFGN$Gene)
# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(genes_cebi_atel, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>% 
  filter(!duplicated(SYMBOL))

genes_cebi_lemu <- unique(cebi_lemu_byFGN$Gene)
write.table(genes_cebi_lemu, file.path(resultsDir,"/functional/genes_cebi_lemu005.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(genes_cebi_lemu, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>% 
  filter(!duplicated(SYMBOL))

genes_lemu_atel <- unique(lemu_atel_byFGN$Gene)
# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(genes_lemu_atel, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

genes_cerco_atel <- unique(cerco_atel_byFGN$Gene)
# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(genes_cerco_atel, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

genes_cerco_cebi <- unique(cerco_cebi_byFGN$Gene)
# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(genes_cerco_cebi, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>%
  filter(!duplicated(SYMBOL))

genes_cerco_lemu <- unique(cerco_lemu_byFGN$Gene)
# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(genes_cerco_lemu, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

genes_shared_all2fam <- unique(genes_shared_all2fam)
# Obtain desired types of IDs from the symbol IDs from 70 genes in common to all contrasts
ann <- bitr(genes_shared_all2fam, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

# Do the same for all genes studied
background_ann <- bitr(unique(background_genes), fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>% 
  filter(!duplicated(SYMBOL))

ego <- enrichGO(gene          = unique(ann$ENTREZID),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "BP",
                pAdjustMethod = "BH",
                minGSSize = 15,
                maxGSSize = 500,
                pvalueCutoff  = 1, 
                qvalueCutoff  = 0.05, 
                readable = TRUE, 
                universe = unique(background_ann$ENTREZID))

dim(ego)
# Eliminate redundant terms, with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
dim(ego)
ego_sign <- as.data.frame(ego)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
ego_sign <- ego_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)

write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_ceat_leche.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_cebi_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_cebi_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_lemu_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_cerco_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_cerco_cebi.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_cerco_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_common_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
ego_data <- ego@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(ego_data[1:15,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

p1 <- goplot(ego)
ggsave(plot = p1, filename = file.path(resultsDir,paste0("functional/DAG_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

# Create dot plot of 20 most significant BP GO terms
barplot(ego, showCategory = 15) 
p2 <- dotplot(ego, showCategory = 15)
ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/dotplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

# Create Emap plot of 15 most significant GO terms
egop <- enrichplot::pairwise_termsim(ego)
p3 <- emapplot(egop, cex.params = list(category_label = 0.7), showCategory = 15)
ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/emapplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

# Cenet plot of 10 most significant BP GO terms
p4 <- cnetplot(ego, showCategory = 10, cex.params = list(category_label = 0.8),cex_label_gene = 0.5)
ggsave(plot = p4, filename = file.path(resultsDir,paste0("functional/cnetplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

ego <- enrichGO(gene          = unique(ann$ENTREZID),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "MF",
                pAdjustMethod = "BH",
                minGSSize = 15,
                maxGSSize = 500,
                pvalueCutoff  = 1, 
                qvalueCutoff  = 0.05, 
                readable = TRUE, 
                universe = unique(background_ann$ENTREZID))

dim(ego)
# Eliminate redundant terms, with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
dim(ego)
ego_sign <- as.data.frame(ego)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
ego_sign <- ego_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)

write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOMF_ceat_leche.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOMF_cebi_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOMF_cebi_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOMF_lemu_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOMF_cerco_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOMF_cerco_cebi.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOMF_cerco_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
ego_data <- ego@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(ego_data[1:15,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

#dim(ego)
p1 <- goplot(ego)
ggsave(plot = p1, filename = file.path(resultsDir,paste0("functional/DAG_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

# Create dot plot of 20 most significant BP GO terms
barplot(ego, showCategory = 10) 
p2 <- dotplot(ego, showCategory = 10)
ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/dotplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

# Create Emap plot of 15 most significant GO terms
egop <- enrichplot::pairwise_termsim(ego)
p3 <- emapplot(egop, cex.params = list(category_label = 0.7), showCategory = 15)
ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/emapplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

# Cenet plot of 10 most significant BP GO terms
p4 <- cnetplot(ego, showCategory = 15, cex.params = list(category_label = 0.8),cex_label_gene = 0.5)
ggsave(plot = p4, filename = file.path(resultsDir,paste0("functional/cnetplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

set.seed(123)
kk <- enrichKEGG(gene         = unique(ann$ENTREZID),
                 organism     = 'hsa',
                 universe = as.character(unique(background_ann$ENTREZID)),
                 pvalueCutoff = 1,
                 qvalueCutoff  = 0.05, 
                 minGSSize     = 15,
                 maxGSSize     = 500)
dim(kk)
kk_sign <- as.data.frame(kk)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_sign <- kk_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)
#write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_ceat_leche.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_cebi_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_cebi_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_lemu_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_cerco_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_cerco_cebi.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_cerco_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_common_genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_data <- kk@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(kk_data[1:4,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

barplot(kk, showCategory = 10) 
dotplot(kk, showCategory = 10)
kkp <- enrichplot::pairwise_termsim(kk)
emapplot(kkp, cex.params = list(category_label = 0.8), showCategory = 4)
cnetplot(kk, showCategory = 10, cex.params = list(category_label = 0.8))

kk <- enrichMKEGG(gene = unique(ann$ENTREZID),
                   universe = unique(background_ann$ENTREZID),
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 0.05,
                   minGSSize     = 15,
                   maxGSSize     = 500)
dim(kk)
# REACTOME ORA
kk <- enrichPathway(gene = unique(ann$ENTREZID),
                    universe = unique(background_ann$ENTREZID), 
                    pvalueCutoff = 1,
                    qvalueCutoff = 0.05,
                    minGSSize     = 15,
                    maxGSSize     = 500,
                    readable = TRUE)
dim(kk)
kk_sign <- as.data.frame(kk)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_sign <- kk_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)
#write.table(ego_sign, file.path(resultsDir,"/functional/REACTOME_ceat_leche.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/REACTOME_cebi_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/REACTOME_cebi_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/REACTOME_lemu_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/REACTOME_cerco_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/REACTOME_cerco_cebi.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/REACTOME_cerco_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

kk <- enrichDO(gene          = unique(ann$ENTREZID),
               ont           = "DO",
                pvalueCutoff  = 1,
                pAdjustMethod = "BH",
                universe      = unique(background_ann$ENTREZID),
                minGSSize     = 15,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
dim(kk)
kk_sign <- as.data.frame(kk)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_sign <- kk_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)
#write.table(ego_sign, file.path(resultsDir,"/functional/DO_ceat_leche.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DO_cebi_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DO_cebi_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DO_lemu_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DO_cerco_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DO_cerco_cebi.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DO_cerco_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_data <- kk@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(kk_data[1:5,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))


kk <- enrichNCG(unique(ann$ENTREZID),
                 pvalueCutoff  = 1,
                 pAdjustMethod = "BH",
                 universe      = unique(background_ann$ENTREZID),
                 minGSSize     = 15,
                 maxGSSize     = 500,
                 qvalueCutoff  = 0.05,
                 readable = TRUE) 
dim(kk)
kk_sign <- as.data.frame(kk)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_sign <- kk_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)
#write.table(ego_sign, file.path(resultsDir,"/functional/DO_ceat_leche.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/NCG_cebi_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/NCG_cebi_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/NCG_lemu_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/NCG_cerco_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/NCG_cerco_cebi.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/NCG_cerco_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_data <- kk@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(kk_data[1:5,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

kk <- enrichDGN(unique(ann$ENTREZID),
                 pAdjustMethod = "BH",
                 universe      = unique(background_ann$ENTREZID),
                minGSSize     = 15,
                maxGSSize     = 500,
                pvalueCutoff  = 1,
                qvalueCutoff  = 0.05,
                readable = TRUE) 
dim(kk)
kk_sign <- as.data.frame(kk)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_sign <- kk_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)

#write.table(ego_sign, file.path(resultsDir,"/functional/DGN_ceat_leche.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DGN_cebi_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DGN_cebi_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DGN_lemu_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DGN_cerco_atel.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DGN_cerco_cebi.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/DGN_cerco_lemu.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_data <- kk@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(kk_data[1:5,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)

kk <- enricher(unique(ann$ENTREZID),
               pAdjustMethod = "BH", 
               universe      = unique(background_ann$ENTREZID),
               minGSSize     = 15,
               maxGSSize     = 500,
               pvalueCutoff  = 1,
               qvalueCutoff  = 0.05,
               TERM2GENE = m_t2g)
dim(kk)
#C7: immunologic signature gene sets
m_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)

kk <- enricher(unique(ann$ENTREZID),
               pAdjustMethod = "BH", 
               universe      = unique(background_ann$ENTREZID),
               minGSSize     = 15,
               maxGSSize     = 500,
               pvalueCutoff  = 1,
               qvalueCutoff  = 0.05,
               TERM2GENE = m_t2g)
dim(kk)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_data <- kk@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(kk_data[1:6,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

ageannots <- read.gmt("data/agein.Hs.symbols.gmt")
kk <- enricher(unique(ann$SYMBOL),
               pAdjustMethod = "BH", 
               universe      = unique(background_ann$SYMBOL),
               minGSSize     = 10,
               maxGSSize     = 500,
               pvalueCutoff  = 1,
               qvalueCutoff  = 0.05,
               TERM2GENE = ageannots)
dim(kk)
kk_sign <- as.data.frame(kk)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_sign <- kk_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)

# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_data <- hpo_results@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(kk_data[1:5,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

