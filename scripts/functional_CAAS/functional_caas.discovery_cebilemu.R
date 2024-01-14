# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(gt)
library(org.Hs.eg.db)
library(clusterProfiler) 
# DOSE supports enrichment analysis of Disease Ontology (DO) (Schriml et al. 2011), Network of Cancer Gene (A. et al. 2016) and Disease Gene Network (DisGeNET) (Janet et al. 2015).
library(DOSE)
library(ReactomePA)
library(msigdbr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(magick)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load data 
cebi_lemu_byFGN.discov <- read_delim(file.path(resultsDir,"/functional/cebi_lemu_discovery.txt"), 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
cebi_lemu.bootstrap <- read_delim(file.path(resultsDir,"caas/bootstrap_results/all.boot"), 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE, col_names = F)

# Background genes
background_genes <- read.table(file.path(resultsDir,"/caas/caastools_LQ_out/background_genes.txt"), header = FALSE, stringsAsFactors = FALSE)$V1

# Name columns in  bootstrap file
colnames(cebi_lemu.bootstrap) <- c("GenePos", "Counts", "Cycles", "adj.pval", "Bootsrap_greater", "Trait")
# keep only significant results in bootstrap file
cebi_lemu.bootstrap <- cebi_lemu.bootstrap %>%
  filter(adj.pval <= 0.05) 

# Make a column as in bootstrap file
cebi_lemu_byFGN.discov <- cebi_lemu_byFGN.discov %>%
  mutate(GenePos = paste(Gene, Position, sep = "@"))

# Perform left join discovery dataframe with the info about their bootstrap pvalue from bootstrap dataframe
longevity_discovery <- left_join(cebi_lemu_byFGN.discov, cebi_lemu.bootstrap[, c("GenePos", "adj.pval", "Counts")],  by = "GenePos")

# keep caas with  significant bootstrap pvalue criteria
longevity_bootstrap.05 <- longevity_discovery %>%
  filter(adj.pval <= 0.05)

# Count the number of CAAS for each unique adj.pval
cebi_lemu_bootstrap_counts <-  longevity_discovery %>% #
  filter(adj.pval <= 0.05) %>%
  group_by(adj.pval) %>%
  summarize(count = n()) %>%               # Count occurrences of each adj.pval
  ungroup() %>%
  arrange(adj.pval) %>%
  mutate(cumulative_count = cumsum(count)) # Calculate the cumulative count

# Calculate the total number of CAAS occurrences for normalization
total_count <- sum(cebi_lemu_bootstrap_counts$count)

# Calculate the probability of observing more than 2927 occurrences
# Assuming each row corresponds to a unique CAAS
p_obs_more_than_2927 <- sum(cebi_lemu_bootstrap_counts$count[cebi_lemu_bootstrap_counts$count > 2927]) / total_count

# Create the density plot
p <- ggplot(cebi_lemu_bootstrap_counts, aes(x = count, y =  ..count../total_count)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 2927, linetype = "dashed", color = "red") +
  labs(x = "CAAS", y = "") +# Remove y-axis title
  theme_classic()  # Use a classic theme

# Print the plot and probability
p
print(paste("P(X > obs) =", p_obs_more_than_2927))
# # Calculate the density
# density_data <- density(cebi_lemu_bootstrap_counts$count)
# 
# # Integrate the tail of the density from obs to infinity (or the maximum observed count)
# p_obs_more_than_2927 <- sum(density_data$y[density_data$x > 2927]) * diff(range(density_data$x))/length(density_data$x)
# 
# # Display the probability
# print(paste("P(X > obs) =", p_obs_more_than_2927))

# Create the histogram with percentages on the y-axis
p <- ggplot(cebi_lemu_bootstrap_counts, aes(x = count, y = (..count..)/total_count )) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 2927, linetype = "dashed", color = "red") +
  scale_y_continuous() +
  theme_classic()  # Use a minimal theme

# Print the plot
p


# # Perform left join nominal significant dataframe with the info about their bootstrap pvalue from bootstrap dataframe
# longevity_nominal <- left_join(longevity_nominal.05, cebi_lemu.bootstrap[, c("GenePos", "adj.pval")],  by = "GenePos")
# longevity_bootstrap.05 <- longevity_nominal %>%
#   filter(adj.pval <= 0.05)

# Create a table like the one in Farre's paper:
# Summarize the Data
summary_data <- longevity_discovery %>%
  filter(Pvalue <= 0.05)  %>%
  group_by(Pattern) %>%
  # Summarize discovery data
  summarize(
    CAAS_discov = n(),
    Genes_discov = n_distinct(Gene),
    CAAS_valid = sum(adj.pval <= 0.05, na.rm = TRUE), # Count instances where adj.pval is <= 0.05
    # Count unique genes where adj.pval is <= 0.05
    Genes_valid = n_distinct(Gene[!is.na(adj.pval) & adj.pval <= 0.05])) %>% 
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
    Pattern = "Pattern",
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

# Save the gt table as a Word document (experimental in `gt`)
gt::gtsave(table_display, filename = file.path(resultsDir, "functional/caas_table.docx"))

# Save unique symbol IDs
genes_cebi_lemu <- unique(longevity_bootstrap.05$Gene)
write.table(genes_cebi_lemu, file = file.path(resultsDir,"genes_cebi_lemut.tsv"), sep = "\t",col.names = FALSE, row.names = FALSE, quote = FALSE)

# Obtain desired types of IDs from the symbol IDs from significant genes
#ann <- bitr(genes_cebi_lemu, fromType = "SYMBOL", toType =  c("SYMBOL","ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db) %>% 
#  filter(!duplicated(SYMBOL))

# Export IDs genes with significant bootstrap aa substitutions for using it in S-LDSC
#write.table(as.vector(ann$ENSEMBL), file = file.path(resultsDir,"signif_caastools_geneset.tsv"), sep = "\t",col.names = FALSE, row.names = FALSE, quote = FALSE)

# Keep nominal significant results
longevity_nominal.05 <- longevity_discovery %>%
  filter(Pvalue <= 0.05)
# Obtain unique gene symbols for discovery caas
genes_cebi_lemu <- unique(longevity_nominal.05$Gene)
# Obtain equivalent IDS
ann <- bitr(genes_cebi_lemu, fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>% 
  filter(!duplicated(SYMBOL))

# Do the same for all genes studied
background_ann <- bitr(unique(background_genes), fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>% 
  filter(!duplicated(SYMBOL))

# Check if there are NAs
na_count <- sum(is.na(background_ann$ENTREZID))
# subset dataframe if any NA values are found
if (na_count > 0) {
  background_ann <- subset(background_anns, !is.na(ENTREZID))
  cat("NA values found and rows containing them have been removed.\n")
} 

# OVER REPRESENTATION ANALYSIS (ORA)
  
#  To determine what a priori defined gene sets  are more represented that we could expect by chance in our subset of significant genes 
# Calculate ORA for biological processes (BP) de GO, with minGSSize and maxGSSize to restringe gene sets size
#minGSSize, is minimal number of genes annotated to an Ontology term  that will be tested (to exclude really small ones because they would reach easily significance). maxGSSize is maximum number of annotated genes that will be tested (to descart really big ones, that will be too general/unspecific)

# Define GO terms to study
category <- c("BP", "MF", "CC")

# Initialize empty list 
egos <- list()


for (i in category) {
  set.seed(123)  #Set random seed, so result doesn't change each time
  print(paste("ORA for GO", i))
  ego <- enrichGO(gene          = unique(ann$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = i,
                 pAdjustMethod = "BH",
                 minGSSize = 15,
                 maxGSSize = 500,
                 pvalueCutoff  = 1, 
                 qvalueCutoff  = 0.05, 
                 readable = TRUE, 
                 universe = unique(background_ann$ENTREZID))
  
  # Save ora object in list
  egos[[i]] <- c(egos[[i]], ego)
  # Save results table
  write.csv(ego@result, file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuGO_", i,".csv")), row.names = FALSE) # /cebi_lemubootstrap/
  
  # Number significant BP
  res <- dim(ego) 
  print(res)
}

for (i in 1:length(egos)) {
  if (dim(as.data.frame(egos[[i]]))[1] > 0) {
  ego_sign <- as.data.frame(egos[[i]])
  # Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
  ego_sign <- ego_sign %>%
    separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
    separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
    mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
    mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom),
           ShortDescription = str_wrap(Description, width = 30)) %>% # Create a new column with shortened descriptions
    select(ID, ShortDescription, Description, ER,qvalue)
  
  write.table(ego_sign, file.path(resultsDir,paste0("functional/cebi_lemubootstrap/ERatio_GO_", names(egos)[i],".csv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  # Visualize results graphically, to better understand main processes affected by these genes:
  if (dim(as.data.frame(egos[[i]]))[1] > 10) {
    # Eliminate redundant terms only for BP, , with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
    ego <- simplify(egos[[i]][[1]], cutoff = 0.7, by = "p.adjust", select_fun = min)
    dim(ego)
    ego_sign <- as.data.frame(ego)
    # Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
    ego_sign <- ego_sign %>%
      separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
      separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
      mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
      mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom),
             ShortDescription = str_wrap(Description, width = 30)) %>% # Create a new column with shortened descriptions
      select(ID, ShortDescription, Description, ER,qvalue)
    # Plot significant results as ER barplot
    p <- ggplot(ego_sign[1:10,], aes(x = reorder(ShortDescription, ER), y = ER, fill = -qvalue)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_continuous(low = "blue", high = "red") +
      labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
      theme_classic() +
      theme(axis.text = element_text(size = 11)) 
    ggsave(plot = p, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuERbarplot_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
    
    # Create DAG of significant BP GO terms
    p1 <- goplot(egos[[i]][[1]])
    ggsave(plot = p1, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuDAG_GO_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
    
    # Create dot plot of 20 most significant BP GO terms
    p2 <- dotplot(egos[[i]][[1]], showCategory = 10, font.size = 7)
    ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemudotplot_GO_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
    
    # Create Emap plot of 20 most significant GO terms
    egop <- enrichplot::pairwise_termsim(egos[[i]][[1]])
    p3 <- emapplot(egop, cex.params = list(category_label = 0.8), showCategory = 20)
    ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuemapplot_GO_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
  
    # Cenet plot of 10 most significant BP GO terms
    p4 <- cnetplot(egos[[i]][[1]], showCategory = 7, cex.params = list(category_label = 0.8))
    ggsave(plot = p4, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemucnetplot_GO_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
  } else {
    # Plot significant results as ER barplot
    p <- ggplot(ego_sign, aes(x = reorder(ShortDescription, ER), y = ER, fill = -qvalue)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_continuous(low = "blue", high = "red") +
      labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
      theme_classic() +
      theme(axis.text = element_text(size = 11)) 
    ggsave(plot = p, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuERbarplot_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
    
    # Create DAG of significant BP GO terms
    p1 <- goplot(egos[[i]][[1]])
    ggsave(plot = p1, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuDAG_GO_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
    
    # Create dot plot of 20 most significant BP GO terms
    p2 <- dotplot(egos[[i]][[1]], showCategory = dim(as.data.frame(egos[[i]]))[1], font.size = 5)
    ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemudotplot_GO_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
    
    # Create Emap plot of 20 most significant GO terms
    egop <- enrichplot::pairwise_termsim(egos[[i]][[1]])
    p3 <- emapplot(egop, cex.params = list(category_label = 0.6), showCategory = dim(as.data.frame(egos[[i]]))[1])
    ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuemapplot_GO_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
    
    # Cenet plot of 10 most significant BP GO terms
    p4 <- cnetplot(egos[[i]][[1]], showCategory = dim(as.data.frame(egos[[i]]))[1], cex.params = list(category_label = 0.6))
    ggsave(plot = p4, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemucnetplot_GO_", names(egos)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
  }
  }
}

# Eliminate redundant terms, with a similarity > than 0.7 and select as representative term from the redundant ones the term with the smallest adjusted pvale
#ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
# s_ego <- clusterProfiler::simplify(egomf) # same results as above
# Number significant BP GOs 
#after <- dim(ego)
#message <- sprintf("Terms significant before simplification: %s. Terms significant after: %s", before[1], after[1])
#print(message)


# prepare figure 4
# Read the images
image1 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemuERbarplot_BP.tiff"))
image2 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemuERbarplot_MF.tiff"))
image3 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemuERbarplot_CC.tiff"))

# Annotate each image with letters A, B, and C respectively
image1 <- ggdraw() + draw_image(image1)
image2 <- ggdraw() + draw_image(image2)
image3 <- ggdraw() + draw_image(image3)

# Combine the annotated images
combined_plot <- image1 + image2 + image3 + 
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = 'A')  # Arrange in one row with 3 columns

combined_plot <- (image1 + image2) / image3  + 
  plot_annotation(tag_levels = 'A')
combined_plot

# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"functional/ORA_results/cebi_lemuERbarplot_GO.png"),  dpi = 300)

# Prepare figure 5
# Read the images
image1 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemuemapplot_GO_BP.tiff"))
image2 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemuemapplot_GO_MF.tiff"))
image3 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemuemapplot_GO_CC.tiff"))

# Annotate each image with letters A, B, and C respectively
image1 <- ggdraw() + draw_image(image1)
image2 <- ggdraw() + draw_image(image2)
image3 <- ggdraw() + draw_image(image3)

# Combine the annotated images
combined_plot <- image1 / (image2 + image3) + 
  plot_layout(heights = c(2, 1.5)) + 
  plot_annotation(tag_levels = 'A')
combined_plot
# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"functional/ORA_results/cebi_lemuemapplot_GO.png"),  dpi = 300)


# KEGG pathway over-representation analysis
set.seed(123)
kk <- enrichKEGG(gene         = unique(ann$ENTREZID),
                 organism     = 'hsa',
                 universe = as.character(unique(background_ann$ENTREZID)),
                 pvalueCutoff = 1,
                 qvalueCutoff  = 0.05,
                 minGSSize     = 15,
                 maxGSSize     = 500)

#KEGG module over-representation analysis
#KEGG Module is a collection of manually defined function units. In some situation, KEGG Modules have a more straightforward interpretation
mkk <- enrichMKEGG(gene = unique(ann$ENTREZID),
                   universe = unique(background_ann$ENTREZID),
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 0.05,
                   minGSSize     = 15,
                   maxGSSize     = 500)
# REACTOME ORA
pkk <- enrichPathway(gene = unique(ann$ENTREZID),
                    universe = unique(background_ann$ENTREZID), 
                    pvalueCutoff = 1,
                    qvalueCutoff = 0.05,
                    minGSSize     = 15,
                    maxGSSize     = 500,
                    readable = TRUE)

# Disease over-representation analysis
edo <- enrichDO(gene          = unique(ann$ENTREZID),
                ont           = "DO",
                pvalueCutoff  = 1,
                pAdjustMethod = "BH",
                universe      = unique(background_ann$ENTREZID),
                minGSSize     = 15,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# Over-representation analysis for the network of cancer gene
# Network of Cancer Gene (NCG) (A. et al. 2016) is a manually curated repository of cancer genes. NCG release 5.0 (Aug. 2015) collects 1,571 cancer genes from 175 published studies. DOSE supports analyzing gene list and determine whether they are enriched in genes known to be mutated in a given cancer type.
ncg <- enrichNCG(unique(ann$ENTREZID),
                 pvalueCutoff  = 1,
                 pAdjustMethod = "BH",
                 universe      = unique(background_ann$ENTREZID),
                 minGSSize     = 15,
                 maxGSSize     = 500,
                 qvalueCutoff  = 0.05,
                 readable = TRUE)  

# Over-representation analysis for the disease gene network
#DisGeNET(Janet et al. 2015) is an integrative and comprehensive resources of gene-disease associations from several public data sources and the literature. It contains gene-disease associations and snp-gene-disease associations.
#The enrichment analysis of disease-gene associations is supported by the enrichDGN function and analysis of snp-gene-disease associations is supported by the enrichDGNv function.
dgn <- enrichDGN(unique(ann$ENTREZID),
                 pAdjustMethod = "BH",
                 universe      = unique(background_ann$ENTREZID),
                 minGSSize     = 15,
                 maxGSSize     = 500,
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 0.05,
                 readable = TRUE) 

# MSigDb analysis
# we retrieve the dataset: C3: motif gene sets, Gene sets representing potential targets of regulation by transcription factors or microRNAs
#These gene sets make it possible to link changes in an expression profiling experiment to a putative cis-regulatory element. The C3 collection is divided into two subcollections: microRNA targets (MIR) and transcription factor targets (TFT).
m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)

em <- enricher(unique(ann$ENTREZID),
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

emc7 <- enricher(unique(ann$ENTREZID),
               pAdjustMethod = "BH", 
               universe      = unique(background_ann$ENTREZID),
               minGSSize     = 15,
               maxGSSize     = 500,
               pvalueCutoff  = 1,
               qvalueCutoff  = 0.05,
               TERM2GENE = m_t2g)
ageannots <- read.gmt("data/agein.Hs.symbols.gmt")
ageingora <- enricher(unique(ann$SYMBOL),
               pAdjustMethod = "BH", 
               universe      = unique(background_ann$SYMBOL),
               minGSSize     = 10,
               maxGSSize     = 500,
               pvalueCutoff  = 1,
               qvalueCutoff  = 0.05,
               TERM2GENE = ageannots)

# Gather all ORA results in a vector w/ their names
ora_objects <- c(kk,mkk,pkk,edo,ncg,dgn,em,emc7,ageingora)
names(ora_objects) <- c("KEGG", "KEGG_module", "REACTOME", "DOSE", "NCG", "DisGeNET", "MSigDbC3","MSigDbC7", "Ageing")

# For every ORA result
for (i in 1:length(ora_objects)) {
  if (dim(ora_objects[[i]])[1] > 0) {
    cat(paste("\nSaving", names(ora_objects)[i],"results table\n"))
    write.csv(ora_objects[[i]]@result, file.path(resultsDir,paste0("functional/ORA_results/cebi_lemu", names(ora_objects)[i],".csv")), row.names = FALSE)
    
    sign <- as.data.frame(ora_objects[[i]])
    # Convert GeneRatio and BgRatio to numeric and calculate Enrichment Ratio
    sign <- sign %>%
      separate(GeneRatio, into = c("GeneNum", "GeneDenom"), sep = "/") %>%
      separate(BgRatio, into = c("BgNum", "BgDenom"), sep = "/") %>%
      mutate(across(c(GeneNum, GeneDenom, BgNum, BgDenom), as.numeric)) %>%
      mutate(ER = (GeneNum / GeneDenom) / (BgNum / BgDenom)) %>%
      mutate(ShortDescription = str_wrap(Description, width = 30)) %>%
      select(ID, ShortDescription, Description, ER, qvalue)
      cat(paste("Saving", names(ora_objects)[i]," Enrichment Ration results table\n"))
      write.table(sign, file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuERatio", names(ora_objects)[i],".tsv")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
      
      cat(paste("Creating plots for", names(ora_objects)[i],"results\n"))
      if (dim(ora_objects[[i]])[1] > 10) {
        # plot significant results as ER barplot
        p <- ggplot(sign[1:10,], aes(x = reorder(ShortDescription, ER), y = ER, fill = -qvalue)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          scale_fill_continuous(low = "blue", high = "red") +
          labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
          theme_classic() +
          theme(axis.text = element_text(size = 11))
        ggsave(plot = p, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuERbarplot_", names(ora_objects)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
        
        # Create dot plot of 10 most significant  terms
        p1 <- dotplot(ora_objects[[i]], showCategory = 10, font.size = 8)
        ggsave(plot = p1, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemudotplot_", names(ora_objects)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
      
        # Create Emap plot of 10 most significant terms
        ekk <- enrichplot::pairwise_termsim(ora_objects[[i]])
        p2 <- emapplot(ekk, showCategory = 10)
        ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuemapplot_", names(ora_objects)[i],".png")), dpi = 300,units = "in",width = 8, height = 8)
      
        # Cenet plot of 10 most significant terms
        p3 <- cnetplot(ora_objects[[i]], showCategory = 7, cex.params = list(category_label = 0.8))
        ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemucnetplot_", names(ora_objects)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
        } else {
          # plot significant results as ER barplot
          p <- ggplot(sign, aes(x = reorder(ShortDescription, ER), y = ER, fill = -qvalue)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          scale_fill_continuous(low = "blue", high = "red") +
          labs(x = "", y = "Enrichment Ratio", fill = "FDR") +
          theme_classic() +
          theme(axis.text = element_text(size = 11))
          ggsave(plot = p, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuERbarplot_", names(ora_objects)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
      
          # Create dot plot of 10 most significant  terms
          p1 <- dotplot(ora_objects[[i]], showCategory = dim(ora_objects[[i]])[1], font.size = 8)
          ggsave(plot = p1, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemudotplot_", names(ora_objects)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
          
          # Create Emap plot of 10 most significant terms
          ekk <- enrichplot::pairwise_termsim(ora_objects[[i]])
          p2 <- emapplot(ekk, showCategory = dim(ora_objects[[i]])[1])
          ggsave(plot = p2, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemuemapplot_", names(ora_objects)[i],".png")), dpi = 300,units = "in",width = 8, height = 8)
          
          # Cenet plot of 10 most significant terms
          p3 <- cnetplot(ora_objects[[i]], showCategory = 10, cex.params = list(category_label = 0.6))
          ggsave(plot = p3, filename = file.path(resultsDir,paste0("functional/ORA_results/cebi_lemucnetplot_", names(ora_objects)[i],".png")), dpi = 300, units = "in",width = 6, height = 6)
        }
      } 
  else {
    cat(paste("\nNot significant results found for", names(ora_objects)[i] ,"ORA\n"))
  }
}

# Prepare figure 7
# Read the images
image1 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemudotplot_DOSE.png"))
image2 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemudotplot_DisGeNET.png"))
image3 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemuemapplot_DOSE.png"))
image4 <- magick::image_read(file.path(resultsDir,"functional/ORA_results/cebi_lemuemapplot_DisGeNET.png"))

# Annotate each image with letters A, B, and C respectively
image1 <- ggdraw() + draw_image(image1)
image2 <- ggdraw() + draw_image(image2)
image3 <- ggdraw() + draw_image(image3)
image4 <- ggdraw() + draw_image(image4)

# Combine the annotated images
combined_plot <- (image1 | image2) / (image3 | image4) + 
  plot_annotation(tag_levels = 'A')
combined_plot
# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"functional/ORA_results/cebi_lemudotplotemapplot_DO.png"),  dpi = 300)

  
# Visualize S-LDSC
###################
# Transform backgroun gene symbols to Ensembl IDs to use them in the S-LDSC scripts
background_ids <- bitr(unique(background_genes), fromType = "SYMBOL", toType =  c("SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db) %>% 
  filter(!duplicated(SYMBOL))
# Export them
write.table(background_ids$ENSEMBL, file.path(dataDir,"allgenes.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# Read the data
data1 <- read.table(file.path(resultsDir,"sldsc_out/longev_caastools_geneset_baselineLD.results"), header = TRUE)
# Save results
write.table(data, file.path(resultsDir,"sldsc_out/cebi_lemu_sldsc"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

data <- data1 %>%
  filter(str_starts(Category, regex("^(b|L|Conserved)")) & !str_starts(Category, "Conserved_Lindblad"))

# Create Barplot
barp <- data %>% 
  arrange(desc(Prop._h2)) %>% 
  mutate(Category = str_remove(Category, "_0")) %>% 
  mutate(Category = factor(Category, level = Category)) %>% 
  pivot_longer(.,cols = c(Prop._SNPs,Prop._h2),names_to = "Proportion") %>% 
  ggplot(., aes(x = Category, y = value)) + 
  geom_bar(aes(fill = Proportion),stat = "identity",position = "dodge") +
  theme_classic() +
  ggtitle("Proportion of h2 explained and SNPs used by each category") +
  ylab("proportion") + xlab("Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggsave(plot = barp, filename = file.path(resultsDir,"sldsc_out/sldsc_barplot.tiff"), dpi = 300, units = "in",width = 6, height = 6)

# enrichment plot
# The dotted line shows the bonferonni significance at α cut off of 0.05.
enrchp <- data %>% 
  mutate(Category = str_remove(Category, "_0")) %>% 
  ggplot(., aes(x = Category, y = -log10(Enrichment_p))) +
  geom_hline(yintercept = -log10(0.05/nrow(data1)),linetype = 2) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() + 
  ggtitle("Enrichment of different categories") +
  ylab("-log10(p)") + xlab("Category") +
  coord_flip()

ggsave(plot = enrchp, filename = file.path(resultsDir,"sldsc_out/sldsc_enrchplot.tiff"), dpi = 300, units = "in",width = 6, height = 6)

# Clean h2 results

d <- tibble(filename = data1) %>% mutate(file_contents = map(filename,read_tsv)) %>%
  mutate(makenames = gsub(".results","",basename(filename)),
         Trait = sub("_.*","",makenames),
         File = basename(filename)) %>% 
  select(-filename,-makenames) %>%
  unnest(file_contents) %>%
  mutate(Coefficient_P = pnorm(`Coefficient_z-score`,lower.tail = FALSE)) %>%
  select(File,Trait,Category,Prop._SNPs:`Coefficient_z-score`,Coefficient_P) %>%
  mutate(Category = gsub(".bedL2_0|L2_1|L2_0","",Category)) %>%
  arrange(Trait,Enrichment_p)
  
write_tsv(d,path = "all.results.txt")

sessionInfo()
