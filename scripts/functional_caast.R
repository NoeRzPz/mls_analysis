# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler, quietly = TRUE) 

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load data 
# discovery data frame 
longevity_farre.discov <- read_delim(file.path(dataDir,"longevity.farre2021.nofilter.tab"), 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
# guided bootstrap data frame
bootstrap_df <- read_delim(file.path(dataDir,"guided.tab"), delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE, col_names = FALSE)
  

# Keep nominal significant results
longevity_nominal.05 <- longevity_farre.discov %>%
  filter(Pvalue <= 0.05) %>%
  mutate(id2 = paste(Gene, Position, sep = "@"))

# Perform left join nominal significant dataframe with the info about their bootstrap pvalue from bootstrap dataframe
longevity__nominal <- left_join(longevity_nominal.05, bootstrap_df[, c("X1", "X4")],  by = c("id2" = "X1"))

# Keep bootstrap significant results
longevity_nominal.05 <- longevity__nominal %>%
  rename(adj.pval = X4) %>%
  filter(adj.pval <= 0.05)

# Save unique RefSeq IDs
refseq_id <- unique(longevity_nominal.05$Gene)

# Obtain desired types of IDs from the RefSeq IDs from significant genes
ann <- bitr(refseq_id, fromType = "REFSEQ", toType =  c("REFSEQ","ENSEMBL","ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)

# Do the same for all genes studied
background_genes <- bitr(unique(longevity_farre.discov$Gene), fromType = "REFSEQ", toType =  c("REFSEQ","ENSEMBL","ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)

# Check if there are NAs
na_count <- sum(is.na(background_genes$ENTREZID))
# subset dataframe if any NA values are found
if (na_count > 0) {
  background_genes <- subset(background_genes, !is.na(ENTREZID))
  cat("NA values found and rows containing them have been removed.\n")
} 

# OVER REPRESENTATION ANALYSIS (ORA)
  
#  To determine what a priori defined gene sets  are more represented that we could expect by chance in our subset of significant genes 

# We provide ENSEMBL IDs of our significant genes.  Specify option readable TRUE to translate ENSEMBL  IDs to gene symbols

# Calculate ORA for biological processes (BP) de GO, with minGSSize and maxGSSize to restringe gene sets size
#minGSSize, is minimal number of genes annotated to an Ontology term  that will be tested (to exclude really small ones because they would reach easily significance). maxGSSize is maximum number of annotated genes that will be tested (to descart really big ones, that will be too general/unspecific)

set.seed(123)  #Set random seed, so result doesn't change each time

ego1 <- enrichGO(gene          = unique(ann$ENSEMBL),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 minGSSize = 15,
                 maxGSSize = 500,
                 pvalueCutoff  = 0.05, #0.01
                 qvalueCutoff  = 0.05,
                 readable = TRUE, 
                 universe = background_genes$ENSEMBL)
# Number significant BP
dim(ego1) 

# `simplify` method eliminates redundant terms, with a similarity > than 0.7. 'select_fun' selects as representative term from the redundant ones  the term with the smallest adjusted pvalor. 
ego1 <- simplify(ego1 , cutoff = 0.7, by = "p.adjust", select_fun = min)

# Number significant BP GOs 
dim(ego1)

# Capture significant BP 
sig.BP <- head(ego1, nrow(ego1)) 
head(sig.BP)  

# Save results table
write.csv(ego1@result, file.path(resultsDir,"functional/GO_BP.csv"), row.names = FALSE)

# Viaulize results graphically, to better understand main processes affected by these genes:

# Create DAG of significant BP GO terms
p1 <- goplot(ego1)
ggsave(plot = p1, filename = file.path(resultsDir,"functional/DAG_GO_BP.tiff"), dpi = 300, units = "in",width = 6, height = 6)

# Create dot plot of 20 most significant BP GO terms
p2 <- dotplot(ego1, showCategory = 20, font.size = 5)
ggsave(plot = p2, filename = file.path(resultsDir,"functional/dotplot_GO_BP.tiff"), dpi = 300, units = "in",width = 6, height = 6)

# Create Emap plot of 20 most significant GO terms
ego1p <- enrichplot::pairwise_termsim(ego1)
p3 <- emapplot(ego1p, cex.params = list(category_label = 0.5), showCategory = 20)
ggsave(plot = p3, filename = file.path(resultsDir,"functional/emapplot_GO_BP.tiff"), dpi = 300, units = "in",width = 6, height = 6)

# Cenet plot of 10 most significant BP GO terms
p4 <- cnetplot(ego1, showCategory = 10, cex.params = list(category_label = 0.6))
ggsave(plot = p4, filename = file.path(resultsDir,"functional/cnetplot_GO_BP.tiff"), dpi = 300, units = "in",width = 6, height = 6)

# KEGG pathway over-representation analysis
set.seed(123)

kk <- enrichKEGG(gene         = ann$ENTREZID,
                 organism     = 'hsa',
                 universe = background_genes$ENTREZID,
                 pvalueCutoff = 0.05,
                 qvalueCutoff  = 0.05)

# Save KEGG results table
write.csv(kk@result, file.path(resultsDir,"fuctional/kegg.csv"), row.names = FALSE)

# enrichment map
ekk <- enrichplot::pairwise_termsim(kk)
emapplot(ekk, showCategory = 10)
ggsave(plot = last_plot(), filename = file.path(resultsDir,"functional/emapplot_kegg.tiff"), dpi = 300,units = "in",width = 8, height = 8)

#KEGG module over-representation analysis
#KEGG Module is a collection of manually defined function units. In some situation, KEGG Modules have a more straightforward interpretation
mkk <- enrichMKEGG(gene = ann$ENTREZID,
                   universe = background_genes$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.25)
dim(mkk)

# Save significant results
all.mkk <- head(mkk, nrow(mkk)) 
head(all.mkk) 

#enrichment map
mkk3 <- enrichplot::pairwise_termsim(mkk)
emapplot(mkk3, showCategory = 20)

# Ideas: 
# Study if those genes share a specific regulation eg:if they have binding sites to common transcription factors
# Study if they are regulated by the same miRNAs. 

sessionInfo()
