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

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load data 
# discovery data frame 
caas.discovery <- read_delim(file.path(resultsDir,"/caas/caastools_LQ_out/all.caas_discovery.tsv"), 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#ferre.caas <- read_delim(file.path(resultsDir,"/caas/caastools_LQ_out/all.caas_discovery.tsv"), 
                         #delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# Backgroun genes
background_genes <- read.table(file.path(resultsDir,"/caas/caastools_LQ_out/background_genes.txt"), header = FALSE, stringsAsFactors = FALSE)$V1

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
  summarise(Macaca_fascicularis = sum(Macaca_fascicularis),
            Ateles_geoffroyi = sum(Ateles_geoffroyi),
            Sapajus_apella = sum(Sapajus_apella),
            Eulemur_mongoz = sum(Eulemur_mongoz)) %>%
  arrange(FFGN) 

# Remove the FFGN column for the heatmap as it's not needed

# Convert the dataframe to a matrix, excluding the FFGN column
heatmap_data <- as.matrix(freq_spp[,-1])

# Add the FFGN as row names to the matrix
rownames(heatmap_data) <- freq_spp$FFGN

# Now try generating the heatmap
# Escale by number of foreground species found where CAAS were found
pheatmap(heatmap_data, scale = "row", cluster_rows = FALSE,  cluster_cols = TRUE)

# Create the same plot with background spp
freq_spp <- caas.discovery %>%
  mutate(Alouatta_palliata = ifelse(str_detect(FBG, "Alouatta_palliata"), 1, 0),
        Nasalis_larvatus = ifelse(str_detect(FBG, "Nasalis_larvatus"), 1, 0),
        Saguinus_imperator = ifelse(str_detect(FBG, "Saguinus_imperator"), 1, 0),
        Prolemur_simus = ifelse(str_detect(FBG, "Prolemur_simus"), 1, 0)) %>%
  group_by(FBGN) %>%
  summarise(Alouatta_palliata = sum(Alouatta_palliata),
            Nasalis_larvatus = sum(Nasalis_larvatus),
            Saguinus_imperator = sum(Saguinus_imperator),
            Prolemur_simus = sum(Prolemur_simus)) %>%
  arrange(FBGN) 

# Remove the FBGN column for the heatmap as it's not needed
# Convert the dataframe to a matrix, excluding the FFGN column
heatmap_data <- as.matrix(freq_spp[,-1])

# Add the FFGN as row names to the matrix
rownames(heatmap_data) <- freq_spp$FBGN

# Now try generating the heatmap
# Escale by number of foreground species found where CAAS were found
pheatmap(heatmap_data, scale = "row", cluster_rows = FALSE,  cluster_cols = TRUE)




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

# We select gene symbols for the ORA
genes_byFGN <- unique(df_byFGN$Gene)
# Export gene symbols from this subset
write.table(genes_byFGN, file.path(resultsDir,"/functional/genes_byFGN.tab"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# Obtain desired types of IDs from the RefSeq IDs from significant genes
ann <- bitr(symbol_id, fromType = "SYMBOL", toType =  c("SYMBOL","ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)

# Do the same for all genes studied
background_genes <- bitr(unique(background_genes), fromType = "SYMBOL", toType =  c("SYMBOL","ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene          = unique(ann$ENSEMBL),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 600,
                pvalueCutoff  = 1, 
                qvalueCutoff  = 0.25, 
                readable = TRUE, 
                universe = unique(background_genes$ENSEMBL))

dim(ego)
# Eliminate redundant terms, with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
#ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

#dim(ego)
# plot no significant results
ggplot(ego@result[1:15,], aes(x = reorder(Description, Count), y = Count, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "", fill = "FDR") +
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

ego <- enrichGO(gene          = unique(ann$ENSEMBL),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "MF",
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 600,
                pvalueCutoff  = 1, 
                qvalueCutoff  = 0.25, 
                readable = TRUE, 
                universe = unique(background_genes$ENSEMBL))

dim(ego)
# Eliminate redundant terms, with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
#ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

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
                 universe = as.character(unique(background_genes$ENTREZID)),
                 pvalueCutoff = 1,
                 qvalueCutoff  = 0.25, 
                 minGSSize     = 10,
                 maxGSSize     = 600)

barplot(kk, showCategory = 15) 
dotplot(kk, showCategory = 15)
kkp <- enrichplot::pairwise_termsim(kk)
emapplot(kkp, cex.params = list(category_label = 0.8), showCategory = 15)
cnetplot(kk, showCategory = 15, cex.params = list(category_label = 0.8))

kk <- enrichMKEGG(gene = unique(ann$ENTREZID),
                   universe = unique(background_genes$ENTREZID),
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 0.25,
                   minGSSize     = 10,
                   maxGSSize     = 600)

kk <- enrichDO(gene          = unique(ann$ENTREZID),
                ont           = "DO",
                pvalueCutoff  = 1,
                pAdjustMethod = "BH",
                universe      = unique(background_genes$ENTREZID),
                minGSSize     = 10,
                maxGSSize     = 600,
                qvalueCutoff  = 0.25,
                readable      = TRUE)
ggplot(kk@result[1:10,], aes(x = reorder(Description, Count), y = Count, fill = -p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "blue", high = "red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))

kk <- enrichNCG(unique(ann$ENTREZID),
                 pvalueCutoff  = 1,
                 pAdjustMethod = "BH",
                 universe      = unique(background_genes$ENTREZID),
                 minGSSize     = 10,
                 maxGSSize     = 600,
                 qvalueCutoff  = 0.25,
                 readable = TRUE) 

kk <- enrichDGN(unique(ann$ENTREZID),
                 pAdjustMethod = "BH",
                 universe      = unique(background_genes$ENTREZID),
                minGSSize     = 10,
                maxGSSize     = 600,
                pvalueCutoff  = 1,
                qvalueCutoff  = 0.25,
                readable = TRUE) 
m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)

kk <- enricher(unique(ann$ENTREZID),
               pAdjustMethod = "BH", 
               universe      = unique(background_genes$ENTREZID),
               minGSSize     = 10,
               maxGSSize     = 600,
               pvalueCutoff  = 1,
               qvalueCutoff  = 0.25,
               TERM2GENE = m_t2g)

hpo_results <- enrichHP(gene          = unique(ann$ENSEMBL),
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = i,
                        pAdjustMethod = "BH",
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff  = 0.05, 
                        qvalueCutoff  = 0.05, 
                        readable = TRUE, 
                        universe = unique(background_genes$ENSEMBL))







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



