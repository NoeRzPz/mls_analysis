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
library(ggplot2)
library(pheatmap)
library(ape)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load data 
# discovery data frame 
caas.discovery <- read_delim(file.path(resultsDir,"/caas/caastools_LQ_out/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Backgroun genes
background_genes <- read.table(file.path(resultsDir,"/caas/caastools_LQ_out/background_genes.txt"), header = FALSE, stringsAsFactors = FALSE)$V1

# spp families
spp_fam <- read_delim(file.path(dataDir,"/caas_tool_inputs/my_sp2fam_210727.tab"), 
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# CAAS distributions by Number of species found in the FG group (it excludes those ones having an indel)
freq_table <- caas.discovery %>%
  group_by(FFGN) %>%
  summarise(CAAS_N = n()) %>%
  arrange(FFGN)#,
            #Genes = paste(unique(Gene), collapse = ",")

# Export table
write.table(freq_table , file.path(resultsDir,"/functional/CAAS_byFFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
# Create new columns to register if a species is present or not in a CAAS
caas.discovery <- caas.discovery %>%
  mutate(Macaca_fascicularis = ifelse(str_detect(FFG, "Macaca_fascicularis"), 1, 0),
         Ateles_geoffroyi = ifelse(str_detect(FFG, "Ateles_geoffroyi"), 1, 0),
         Sapajus_apella = ifelse(str_detect(FFG, "Sapajus_apella"), 1, 0),
         Eulemur_mongoz = ifelse(str_detect(FFG, "Eulemur_mongoz"), 1, 0))

# Explore what are the species contributing more to each case
freq_spp <- caas.discovery %>%
  group_by(FFGN) %>%
  summarise(`Macaca fascicularis` = sum(Macaca_fascicularis),
            `Ateles geoffroyi` = sum(Ateles_geoffroyi),
            `Sapajus apella` = sum(Sapajus_apella),
            `Eulemur mongoz` = sum(Eulemur_mongoz)) %>%
  arrange(FFGN) 

# Remove the FFGN column for the heatmap as it's not needed

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
unique_families <- unique(annotation$Family)

# Create a color vector where each family is assigned an hexadecimal value
family_colors <- setNames(c("#f8766d", "#9590ff", "#a3a500", "#d89000"), unique_families)

# Create a list for annotation_colors
annotation_colors <- list(Family = family_colors)

# Create heatmap. Escale by number of foreground species found where CAAS were found
pheatmap(heatmap_data, scale = "row", cluster_rows = FALSE,  cluster_cols = TRUE,
         annotation_col = annotation,
         annotation_colors = annotation_colors,
         angle_col = 315, # To show column names with an angle
         fontsize_row = 12) # To show rows with bigger font size

# Create the same plot with background spp
freq_spp <- caas.discovery %>%
  mutate(Alouatta_palliata = ifelse(str_detect(FBG, "Alouatta_palliata"), 1, 0),
        Nasalis_larvatus = ifelse(str_detect(FBG, "Nasalis_larvatus"), 1, 0),
        Saguinus_imperator = ifelse(str_detect(FBG, "Saguinus_imperator"), 1, 0),
        Prolemur_simus = ifelse(str_detect(FBG, "Prolemur_simus"), 1, 0)) %>%
  group_by(FBGN) %>%
  summarise(`Alouatta palliata` = sum(Alouatta_palliata),
            `Nasalis larvatus` = sum(Nasalis_larvatus),
            `Saguinus imperator` = sum(Saguinus_imperator),
            `Prolemur simus` = sum(Prolemur_simus)) %>%
  arrange(FBGN) 

# Remove the FBGN column for the heatmap as it's not needed
# Convert the dataframe to a matrix, excluding the FFGN column
heatmap_data <- as.matrix(freq_spp[,-1])

# Add the FFGN as row names to the matrix
rownames(heatmap_data) <- freq_spp$FBGN

# Ensure the species names match the column names of your heatmap data
heatmap_data_colnames <- colnames(heatmap_data)

# Create the annotation dataframe with family info
annotation <- spp_fam %>%
  mutate(Family = str_replace(Family, "\\w+_", ""),
         Species = str_replace(Species, "_", " ")) %>%
  filter(Species %in% heatmap_data_colnames) %>%
  column_to_rownames(var = "Species")

annotation[rownames(annotation) == "Saguinus imperator",] <- "Cebidae"
# Assuming 'unique_families' is a vector of all unique family names from your 'species_families' data
unique_families <- unique(annotation$Family)

# Create a color vector where each family is assigned an hexadecimal value
family_colors <- setNames(c("#f8766d", "#a3a500", "#9590ff", "#d89000"), unique_families)

# Create a list for annotation_colors
annotation_colors <- list(Family = family_colors)

# Create heatmap. Escale by number of foreground species found where CAAS were found
pheatmap(heatmap_data, scale = "row", cluster_rows = FALSE,  cluster_cols = TRUE,
         annotation_col = annotation,
         annotation_colors = annotation_colors,
         angle_col = 315, # To show column names with an angle
         fontsize_row = 12) # To show rows with bigger font size



# Plot histogram of hipergeometric test pvalues
pva_histo <- ggplot(data = caas.discovery, aes(x = Pvalue)) + 
  geom_histogram(binwidth = 0.02, fill = "blue", color = "black", alpha = 0.7) + 
  theme_classic() +
  labs(title = "Histogram of P-values", x = "P-value", y = "Frequency")

# save plot
ggsave(plot = pva_histo, filename = file.path(resultsDir,"functional/pvals_histo.tiff"), dpi = 300)

# Scatter plot
fg_bg <- caas.discovery %>% 
  ggplot(aes(x = FBGN, y = FFGN, color = FFGN)) +
  geom_jitter() +
  labs( x = 'FBGN',
        y = 'FFGN') +
  theme_classic() + 
  theme(
    legend.position = "bottom",
    text = element_text(size = 16))

# save plot
ggsave(plot = fg_bg, filename = file.path(resultsDir,"functional/fg_bgscatter.tiff"), dpi = 300)

# CAAS distributions by Number of species found in the FG group (it excludes those ones having an indel)
# Count occurrences
freq_table <-  caas.discovery %>%
  count(FBGN, FFGN) %>%
  spread(key = FFGN, value = n, fill = 0)

# Formatting the table using the gt package
formatted_table <- freq_table %>%
  gt() %>%
  tab_header(title = "FFGN") %>%
  cols_label(
    FBGN = "FBGN")

formatted_table

# Export table, buscar como era con gt
write.table(formatted_table, file.path(resultsDir,"/functional/CAAS_byFFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

freq_table <-  caas.discovery %>%
  filter(Pvalue <= 0.05) %>%
  count(FBGN, FFGN) %>%
  spread(key = FFGN, value = n, fill = 0)

# Formatting the table using the gt package
formatted_table <- freq_table %>%
  gt() %>%
  tab_header(title = "FFGN") %>%
  cols_label( FBGN = "FBGN")

# Export table
write.table(formatted_table, file.path(resultsDir,paste0("/functional/4famCAAS_byFFGN005.txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Count occurrences miss
freq_table <-  caas.discovery %>%
  count(MBG, MFG) %>%
  spread(key = MFG, value = n, fill = 0)

# Formatting the table using the gt package
formatted_table <- freq_table %>%
  gt() %>%
  tab_header(title = "MFG") %>%
  cols_label(
    MBG = "MBG")

formatted_table

# Export table, buscar como era con gt
write.table(formatted_table, file.path(resultsDir,"/functional/CAAS_byFFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Count occurrences gap
freq_table <-  caas.discovery %>%
  count(GBG, GFG) %>%
  spread(key = GFG, value = n, fill = 0)

# Formatting the table using the gt package
formatted_table <- freq_table %>%
  gt() %>%
  tab_header(title = "GFG") %>%
  cols_label(
    GBG = "GBG")

formatted_table

# Export table, buscar como era con gt
write.table(formatted_table, file.path(resultsDir,"/functional/CAAS_byFFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


# Scatter plot
gapsfg_bg <- caas.discovery %>% 
  ggplot(aes(x = factor(GBG), y = factor(GFG), color = factor(GFG))) +
  geom_jitter() +
  scale_color_manual(values = c("0" = "darkblue", "1" = "#00BFFF")) +  # Adjust colors here
  labs( x = 'GBG',
        y = 'GFG') +
  theme_classic() + 
  theme(
    legend.position = "bottom",
    text = element_text(size = 16))

# save plot
ggsave(plot = gapsfg_bg, filename = file.path(resultsDir,"functional/gapsfg_bgscatter.tiff"), dpi = 300)

# Scatter plot
missfg_bg <- caas.discovery %>% 
  ggplot(aes(x = MBG, y = MFG, color = MFG)) +
  geom_jitter() +
  labs( x = 'MBG',
        y = 'MFG') +
  theme_classic() + 
  theme(
    legend.position = "bottom",
    text = element_text(size = 16))

# save plot
ggsave(plot = missfg_bg, filename = file.path(resultsDir,"functional/missfg_bgscatter.tiff"), dpi = 300)


# First we performed analysis w/o filtering by nominal pvalue
# Pattern distribution in nominal significant results
freq_table <- caas.discovery %>%
  group_by(Pattern) %>%
  summarise(Frequency = n(),
            Genes = paste(unique(Gene), collapse = ",")) %>%
  arrange(Pattern)

# CAAS filtered by Number of species found in the FG group (it excludes those ones having an indel) and Number species in BG group
df_byFGN <- caas.discovery %>%
  filter(FFGN >= 3  & FBGN >= 3) 

# We identify where are the gaps
gapsdf_byFGN <- df_byFGN %>%
  filter(GFG != 0 | GBG != 0) 
gapsdf_byFGN

# We remove substitutions with gaps
df_byFGN <- df_byFGN %>%
  filter(GFG == 0 & GBG == 0) 

# Count caas and gaps present
tab_byFGN <- df_byFGN %>%
  summarise(CAAS_N = length(Substitution), 
            GFG_N = sum(GFG),
            GBG_N = sum(GBG),
            Pval = sum(Pvalue <= 0.05))
# Export table
write.table(tab_byFGN , file.path(resultsDir,"/functional/CAAS_3FFGN2FBGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

df_byFGN <- df_byFGN %>%
  filter(Pvalue <= 0.05)
write.table(df_byFGN , file.path(resultsDir,"/functional/df_byFGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# We explore the distribution of CAAS by gene, to detect if there are genes that accumulate more CAAS
caas_bygene <- df_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            Position = paste(Position, collapse = ","),
            Substitution = paste(Substitution, collapse = ","))
# Create a one-dimensional table
caasbygene_distribution <- table(caas_bygene$CAAS_N)

# Convert the table to a data frame for plotting
df_caasbygene <- as.data.frame(caasbygene_distribution)

# Calculate the percentages for each category
df_caasbygene <- df_caasbygene %>%
  mutate(Percentage = Freq / sum(Freq) * 100)

# Plot with ggplot2 using the percentages
ggplot(data = df_caasbygene, aes(x = Var1, y = Percentage)) + 
  geom_bar(stat = "identity") +
  labs(x = "Number CAAS/Gene", y = "Percentage") +
  theme_classic()
  
# We select gene symbols for the ORA
genes_byFGN <- unique(df_byFGN$Gene)
# Export gene symbols from this subset
write.table(genes_byFGN, file.path(resultsDir,"/functional/genes_byFGN.tab"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)



# Obtain desired types of IDs from the RefSeq IDs from significant genes
ann <- bitr(genes_byFGN, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

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
ego_sign <- as.data.frame(ego)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
ego_sign <- ego_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOBP_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
ego_data <- ego@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# Eliminate redundant terms, with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
#ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

#dim(ego)
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
p4 <- cnetplot(ego, showCategory = 15, cex.params = list(category_label = 0.8))
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
ego_sign <- as.data.frame(ego)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
ego_sign <- ego_sign %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
  select(ID,Description, ER,p.adjust)
write.table(ego_sign, file.path(resultsDir,"/functional/ORA_results/GOMF_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
# Eliminate redundant terms, with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
#ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
ego_data <- ego@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(ego_data[1:15,], aes(x = reorder(Description, ER), y = ER, fill = p.adjust)) +
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
barplot(ego, showCategory = 15) 
p2 <- dotplot(ego, showCategory = 15)
ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/dotplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

# Create Emap plot of 15 most significant GO terms
egop <- enrichplot::pairwise_termsim(ego)
p3 <- emapplot(egop, cex.params = list(category_label = 0.7), showCategory = 15)
ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/emapplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)

# Cenet plot of 10 most significant BP GO terms
p4 <- cnetplot(ego, showCategory = 15, cex.params = list(category_label = 0.8))
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
write.table(kk_sign, file.path(resultsDir,"/functional/ORA_results/kegg_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
# Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
kk_data <- kk@result %>%
  separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
  separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
  mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
  mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) 

# plot no significant results
ggplot(kk_data[1:3,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

barplot(kk, showCategory = 15) 
dotplot(kk, showCategory = 15)
kkp <- enrichplot::pairwise_termsim(kk)
emapplot(kkp, cex.params = list(category_label = 0.8), showCategory = 15)
cnetplot(kk, showCategory = 15, cex.params = list(category_label = 0.8))

kk <- enrichMKEGG(gene = unique(ann$ENTREZID),
                   universe = unique(background_ann$ENTREZID),
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 0.05,
                   minGSSize     = 15,
                   maxGSSize     = 500)

dim(kk)

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
write.table(kk_sign, file.path(resultsDir,"/functional/REACTOME_results/kegg_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

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
ggplot(kk_data[1:5,], aes(x = reorder(Description, ER), y = ER, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))


fam4.rank <- df_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            pval_mean = mean(Pvalue, na.rm = TRUE)) %>%
  mutate(rank = CAAS_N*(-log(pval_mean)))

# we order genes by CAAS number x -log(pval)
gene_list <- fam4.rank$rank 
# name the vector
names(gene_list) <- as.character(fam4.rank$Gene)
# omit any NA values 
gene_list <- na.omit(gene_list)
# sort the list in decreasing order 
gene_list = sort(gene_list, decreasing = TRUE)

set.seed(123)
gse <- gseGO(geneList = gene_list, 
             ont = "BP", 
             keyType = "SYMBOL",
             minGSSize = 15, 
             maxGSSize = 500,
             pvalueCutoff = 1, 
             verbose = FALSE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH",
             seed = TRUE,
             scoreType = "pos",
             nPermSimple = 1000)
dim(gse)
gse_sign <- as.data.frame(gse) %>%
  filter(qvalue <= 0.05)
dim(gse_sign)

write.table(gse_sign, file.path(resultsDir,"/functional/GSEA_results/GOBP_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

clusterProfiler::ridgeplot(gse, showCategory = 10,core_enrichment = TRUE, label_format = 30)
gseaplot(gse, geneSetID = 1, title = gse$Description[1])

mfgse <- gseGO(geneList = gene_list, 
               ont = "MF", 
               keyType = "SYMBOL",
               minGSSize = 15, 
               maxGSSize = 500,
               pvalueCutoff = 1, 
               verbose = FALSE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH",
               seed = TRUE,
               scoreType = "pos",
               nPermSimple = 1000)
dim(mfgse)
mfgse_sign <- as.data.frame(mfgse) %>%
  filter(qvalue <= 0.05)
dim(mfgse_sign)

write.table(mfgse_sign, file.path(resultsDir,"/functional/GSEA_results/GOMF_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(names(gene_list), fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) 
# Creating a named vector for mapping
mapping <- setNames(ann$ENTREZID, ann$SYMBOL)

# Replace names in gene_list with Entrez IDs
names(gene_list) <- mapping[names(gene_list)]

# Handling any NA values that may have resulted from unmapped symbols
# Remove NA values (unmapped genes)
gene_list <- gene_list[!is.na(names(gene_list))]

kk_gse <- gseKEGG(geneList     = gene_list,
                  organism     = "hsa",
                  minGSSize    = 15, 
                  maxGSSize = 500,
                  pvalueCutoff = 1,
                  nPermSimple = 1000, 
                  pAdjustMethod = "BH",
                  verbose      = FALSE,
                  seed = TRUE,
                  scoreType = "pos")
dim(kk_gse)
kk_gse <- setReadable(kk_gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kkgse_sign <- as.data.frame(kk_gse) %>%
  filter(qvalue <= 0.05)
dim(kkgse_sign)

write.table(kkgse_sign, file.path(resultsDir,"/functional/GSEA_results/KEGG_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

fgsea_react <- gsePathway(geneList = gene_list, 
                          organism = "human",
                          minGSSize = 10, 
                          maxGSSize = 500, 
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          eps = 0, 
                          nPermSimple = 1000, 
                          seed = TRUE,
                          scoreType = "pos")
fgsea_react <- setReadable(fgsea_react, 'org.Hs.eg.db')
dim(fgsea_react)
react_sign <- as.data.frame(fgsea_react) %>%
  filter(qvalue <= 0.05)
dim(react_sign)

write.table(react_sign, file.path(resultsDir,"/functional/GSEA_results/react_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

do_gsea <- gseDO(geneList = gene_list, 
                 minGSSize = 15, 
                 maxGSSize = 500, 
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 nPermSimple = 1000, 
                 seed = TRUE,
                 scoreType = "pos",
                 verbose = FALSE)
do_gsea <- setReadable(do_gsea, 'org.Hs.eg.db')
dim(do_gsea)
do_sign <- as.data.frame(do_gsea) %>%
  filter(qvalue <= 0.05)
dim(do_sign)

write.table(do_sign, file.path(resultsDir,"/functional/GSEA_results/do_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

dgn_gsea <- gseDGN(geneList = gene_list,
                   minGSSize = 15, 
                   maxGSSize = 500, 
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   nPermSimple = 1000, 
                   seed = TRUE,
                   scoreType = "pos",
                   verbose = FALSE)
dgn_gsea <- setReadable(dgn_gsea, 'org.Hs.eg.db')
dim(dgn_gsea)
dgn_sign <- as.data.frame(dgn_gsea) %>%
  filter(qvalue <= 0.05)
dim(dgn_sign)

write.table(dgn_sign, file.path(resultsDir,"/functional/GSEA_results/dgn_4fam.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


hpo_results <- enrichHP(gene          = unique(ann$ENSEMBL),
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = i,
                        pAdjustMethod = "BH",
                        minGSSize = 15,
                        maxGSSize = 500,
                        pvalueCutoff  = 1, 
                        qvalueCutoff  = 0.05, 
                        readable = TRUE, 
                        universe = unique(background_genes$ENSEMBL))

dim(hpo_results)
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

# ferre validated by RRPP 
ferre.discovery <- read_delim(file.path(dataDir,"/ferre_validatedRRPP.txt"),col_names = FALSE ,
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)$X1

ferre.discovery <- unique(ferre.discovery)


# CAAS filtered by Missing species in the FG group and Missing species in BG group
df_byMG <- caas.discovery %>%
  filter(MFG <= 1  & MBG <= 2) 

df_byMG <- df_byMG %>%
  filter(Pvalue <= 0.05)


genes_byMG <- unique(df_byMG$Gene)
# Export gene symbols from this subset
write.table(genes_byMG, file.path(resultsDir,"/functional/genes_byMG.tab"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

tab_byMG <- df_byMG %>%
  summarise(CAAS_N = length(Substitution), 
            GFG_N = sum(GFG),
            GBG_N = sum(GBG))

# Export table
write.table(tab_byMG, file.path(resultsDir,"/functional/CAAS_1MFG2MBG.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# CAAS distributions by Number of species found in the FG group with 3 or more Background species
freq_table <- caas.discovery %>%
  filter(FBGN >= 3) %>%
  group_by(FFGN) %>%
  summarise(Frequency = n(),
            Genes = paste(unique(Gene), collapse = ",")) %>%
  arrange(FFGN)
# Export table
write.table(freq_table , file.path(resultsDir,"/functional/CAAS_byFFGN3FBGN.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


# Keep nominal significant CAAS only
longevity_nominal.05 <- caas.discovery %>%
  filter(Pvalue <= 0.05) #%>%
#  mutate(id2 = paste(Gene, Position, sep = "@"))

# Pattern distribution in nominal significant results
table(longevity_nominal.05$Pattern)
freq_table <- longevity_nominal.05 %>%
  group_by(Pattern) %>%
  summarise(Frequency = n(),
            Genes = paste(unique(Gene), collapse = ",")) %>%
  arrange(Pattern)


# Export filter file
write.table(filtered_caas.discovery, file.path(resultsDir,"/caas/filtered_caas.discovery.tab"), sep = "\t", row.names = FALSE, quote = FALSE)
