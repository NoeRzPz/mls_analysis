# Clean environment
rm(list = ls())

# Load packages
library(readr, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(gt, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)
library(clusterProfiler, quietly = TRUE) 
# DOSE supports enrichment analysis of Disease Ontology (DO) (Schriml et al. 2011), Network of Cancer Gene (A. et al. 2016) and Disease Gene Network (DisGeNET) (Janet et al. 2015).
library(DOSE)
library(msigdbr)
library(ggplot2)

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
longevity_bootstrap.05 <- longevity__nominal %>%
  rename(adj.pval = X4) %>%
  filter(adj.pval <= 0.05)

# Create a table like the one in Farre's paper:
# Summarize the Data
summary_data <- longevity__nominal %>%
  group_by(Scenario) %>%
  # Summarize discovery data
  summarize(
    CAAS_discov = n(),
    Genes_discov = n_distinct(Gene),
    CAAS_valid = sum(X4 <= 0.05, na.rm = TRUE), # Count instances where X4 is <= 0.05
    # Count unique genes where X4 is <= 0.05
    Genes_valid = n_distinct(Gene[X4 <= 0.05])) %>% 
  # Calculate the percentage for validated CAAS and Genes
  mutate(CAAS_valid_perc = CAAS_valid / CAAS_discov * 100,
         Genes_valid_perc = Genes_valid / Genes_discov * 100,
         CAAS_valid = paste0(CAAS_valid, " (", sprintf("%.1f", CAAS_valid_perc), "%)"),
         Genes_valid = paste0(Genes_valid, " (", sprintf("%.1f", Genes_valid_perc), "%)")) %>%
  select(-CAAS_valid_perc, -Genes_valid_perc)

# Formatting the table using gt
table_display <- summary_data %>%
  gt() %>%
  # Set table title
  tab_header(title = "Table. Lists of Discovered and Validated CAAS and Genes") %>%
  # Set column labels
  cols_label(
    Scenario = "Scenario",
    CAAS_discov = "CAAS",
    Genes_discov = "Genes",
    CAAS_valid = "CAAS",
    Genes_valid = "Genes") %>%
  # Group columns under 'Discovered' and 'Validated'
  tab_spanner(
    label = "Discovered",
    columns = c("CAAS_discov", "Genes_discov")) %>%
  tab_spanner(
    label = "Validated",
    columns = c("CAAS_valid", "Genes_valid")) %>%
  # Add calculated numbers and percentages to the table
  fmt_number(
    columns = c("CAAS_discov", "Genes_discov", "CAAS_valid", "Genes_valid"),
    decimals = 0) %>%
  # Add the footnote
  tab_footnote(
    footnote = "NOTE.-Numbers in parentheses represent the percentage of bootstraped validated positions.",
    locations = cells_body(columns = c("CAAS_valid", "Genes_valid")))

# Display the table
table_display

# Save unique RefSeq IDs
refseq_id <- unique(longevity_bootstrap.05$Gene)

# Obtain desired types of IDs from the RefSeq IDs from significant genes
ann <- bitr(refseq_id, fromType = "REFSEQ", toType =  c("REFSEQ","ENSEMBL","ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)

# Export IDs genes with significant bootstrap aa substitutions for using it in S-LDSC
write.table(as.vector(ann$ENSEMBL), file = file.path(resultsDir,"signif_caastools_geneset.tsv"), sep = "\t",col.names = FALSE, row.names = FALSE, quote = FALSE)

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

# ORA 

# Define GO terms to study
category <- c("BP", "MF")

# Initialize empty list (pensar si lo necesito, sino borrar)
egos <- list()

for (i in category) {
  set.seed(123)  #Set random seed, so result doesn't change each time
  print(paste("ORA for GO", i))
  ego <- enrichGO(gene          = unique(ann$ENSEMBL),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = i,
                 pAdjustMethod = "BH",
                 minGSSize = 15,
                 maxGSSize = 500,
                 pvalueCutoff  = 0.05, #0.01
                 qvalueCutoff  = 0.05,
                 readable = TRUE, 
                 universe = unique(background_genes$ENSEMBL))
  
  # Save ora object in list
  egos[[i]] <- c(egos[[i]], ego)
  
  # Number significant BP
  before <- dim(ego) 
  
  # Eliminate redundant terms, with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
  ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  # s_ego <- clusterProfiler::simplify(egomf) # same results as above
  # Number significant BP GOs 
  after <- dim(ego)
   
  message <- sprintf("Terms significant before simplification: %s. Terms significant after: %s", before[1], after[1])
  print(message)
  
  # Save results table
  write.csv(ego@result, file.path(resultsDir,paste0("functional/GO_", i,".csv")), row.names = FALSE)
  
  # Visualize results graphically, to better understand main processes affected by these genes:
  # Create DAG of significant BP GO terms
  p1 <- goplot(ego)
  ggsave(plot = p1, filename = file.path(resultsDir,paste0("functional/DAG_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)
  
  # Create dot plot of 20 most significant BP GO terms
  p2 <- dotplot(ego, showCategory = 20, font.size = 5)
  ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/dotplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)
  
  # Create Emap plot of 20 most significant GO terms
  egop <- enrichplot::pairwise_termsim(ego)
  p3 <- emapplot(egop, cex.params = list(category_label = 0.5), showCategory = 20)
  ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/emapplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)
  
  # Cenet plot of 10 most significant BP GO terms
  p4 <- cnetplot(ego, showCategory = 10, cex.params = list(category_label = 0.6))
  ggsave(plot = p4, filename = file.path(resultsDir,paste0("functional/cnetplot_GO_", i,".tiff")), dpi = 300, units = "in",width = 6, height = 6)
}


# KEGG pathway over-representation analysis
set.seed(123)
kk <- enrichKEGG(gene         = ann$ENTREZID,
                 organism     = 'hsa',
                 universe = as.character(background_genes$ENTREZID),
                 pvalueCutoff = 0.05)

#KEGG module over-representation analysis
#KEGG Module is a collection of manually defined function units. In some situation, KEGG Modules have a more straightforward interpretation
mkk <- enrichMKEGG(gene = ann$ENTREZID,
                   universe = background_genes$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# Disease over-representation analysis
edo <- enrichDO(gene          = ann$ENTREZID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = background_genes$ENTREZID,
              minGSSize     = 15,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = TRUE)

# Over-representation analysis for the network of cancer gene
# Network of Cancer Gene (NCG) (A. et al. 2016) is a manually curated repository of cancer genes. NCG release 5.0 (Aug. 2015) collects 1,571 cancer genes from 175 published studies. DOSE supports analyzing gene list and determine whether they are enriched in genes known to be mutated in a given cancer type.
ncg <- enrichNCG(ann$ENTREZID,
                 pvalueCutoff  = 0.05,
                 pAdjustMethod = "BH",
                 universe      = background_genes$ENTREZID,
                 minGSSize     = 15,
                 maxGSSize     = 500,
                 qvalueCutoff  = 0.05,
                 readable = TRUE) 

# Over-representation analysis for the disease gene network
#DisGeNET(Janet et al. 2015) is an integrative and comprehensive resources of gene-disease associations from several public data sources and the literature. It contains gene-disease associations and snp-gene-disease associations.
#The enrichment analysis of disease-gene associations is supported by the enrichDGN function and analysis of snp-gene-disease associations is supported by the enrichDGNv function.
dgn <- enrichDGN(ann$ENTREZID,
                 pAdjustMethod = "BH",
                 universe      = background_genes$ENTREZID,
                 minGSSize     = 15,
                 maxGSSize     = 500,
                 readable = TRUE) 

# MSigDb analysis
# we retrieve the dataset: C3: motif gene sets, Gene sets representing potential targets of regulation by transcription factors or microRNAs
#These gene sets make it possible to link changes in an expression profiling experiment to a putative cis-regulatory element. The C3 collection is divided into two subcollections: microRNA targets (MIR) and transcription factor targets (TFT).
m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)

em <- enricher(ann$ENTREZID,
               pAdjustMethod = "BH", 
               universe      = background_genes$ENTREZID,
               minGSSize     = 15,
               maxGSSize     = 500,
               TERM2GENE = m_t2g)

# Gather all ORA results in a vector w/ their names
ora_objects <- c(kk,mkk,edo,ncg,dgn,em)
names(ora_objects) <- c("KEGG", "KEGG_module", "Disease", "NCG", "DisGeNET", "MSigDb")

# For every ORA result
for (i in 1:length(ora_objects)) {
  if (dim(ora_objects[[i]])[1] > 0) {
    cat(paste("Saving", names(ora_objects)[i],"results table\n"))
    write.csv(ora_objects[[i]]@result, file.path(resultsDir,paste0("functional/", names(ora_objects)[i],".csv")), row.names = FALSE)
    
    cat(paste("Creating plots for", names(ora_objects)[i],"results\n"))
    # Create dot plot of 10 most significant kegg terms
    p1 <- dotplot(ora_objects[[i]], showCategory = 10, font.size = 5)
    ggsave(plot = p1, filename = file.path(resultsDir,paste0("functional/dotplot_", names(ora_objects)[i],".tiff")), dpi = 300, units = "in",width = 6, height = 6)
    
    # Create Emap plot of 10 most significant kegg terms
    ekk <- enrichplot::pairwise_termsim(ora_objects[[i]])
    p2 <- emapplot(ekk, showCategory = 10)
    ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/emapplot_", names(ora_objects)[i],".tiff")), dpi = 300,units = "in",width = 8, height = 8)
    
    # Cenet plot of 10 most significant kegg terms
    p3 <- cnetplot(ora_objects[[i]], showCategory = 10, cex.params = list(category_label = 0.6))
    ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/cnetplot_", names(ora_objects)[i],".tiff")), dpi = 300, units = "in",width = 6, height = 6)
  } else {
    cat(paste("Not significant results found for", names(ora_objects)[i] ,"ORA\n"))
  }
}


# Visualize S-LDSC
###################
# Read the data
data <- read.table(file.path(resultsDir,"longev_caastools_geneset_baselineLD.results"), header = TRUE)

# Create Barplot
barp <- data %>% 
  arrange(desc(Prop._h2)) %>% 
  mutate(Category = str_remove(Category, "_0")) %>% 
  mutate(Category = factor(Category, level = Category)) %>% 
  pivot_longer(.,cols = c(Prop._SNPs,Prop._h2),names_to = "Proportion") %>% 
  ggplot(., aes(x = Category, y = value)) + 
  geom_bar(aes(fill = Proportion),stat = "identity",position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4)) +
  ggtitle("Proportion of h2 explained and snps used by each category") +
  ylab("proportion") + xlab("Category")

ggsave(plot = barp, filename = file.path(resultsDir,"functional/sldsc_barplot.tiff"), dpi = 300, units = "in",width = 6, height = 6)

# enrichment plot
# The dotted line shows the bonferonni significance at α cut off of 0.05.
enrchp <- data %>% 
  mutate(Category = str_remove(Category, "_0")) %>% 
  ggplot(., aes(x = Category, y = -log10(Enrichment_p))) +
  geom_hline(yintercept = -log10(0.05/nrow(data)),linetype = 2) +
  geom_bar(stat = "identity",position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 0.8)) + 
  ggtitle("Enrichment of different categories") +
  ylab("-log10(p)") + xlab("Category") +
  coord_flip() 
ggsave(plot = enrchp, filename = file.path(resultsDir,"functional/sldsc_enrchplot.tiff"), dpi = 300, units = "in",width = 6, height = 6)

sessionInfo()
