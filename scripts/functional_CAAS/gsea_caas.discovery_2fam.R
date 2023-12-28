# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler) 
library(DOSE)
library(msigdbr)
library(ReactomePA)
library(cowplot)
library(patchwork)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load data 
# discovery data frame 
all4fam.discovery <- read_delim(file.path(resultsDir,"/caas/caastools_LQ_out/all.caas_discovery.tsv"), 
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


# We filter each discovery dataframe contrast by number of FG and BG species in which a CAAS was detected

all4fam_byFGN <- all4fam.discovery %>%
  filter(FFGN >= 3  & FBGN >= 3) 

cebi_atel_byFGN <- cebi_atel.discovery %>%
  filter(FFGN == 2  & FBGN == 2) 

cebi_lemu_byFGN <- cebi_lemu.discovery %>%
  filter(FFGN == 2  & FBGN == 2)

lemu_atel_byFGN <- lemu_atel.discovery %>%
  filter(FFGN == 2  & FBGN == 2) 

cerco_atel_byFGN <- cerco_atel.discovery %>%
  filter(FFGN == 4  & FBGN == 4) 

cerco_cebi_byFGN <- cerco_cebi.discovery %>%
  filter(FFGN == 4  & FBGN == 4) 

cerco_lemu_byFGN <- cerco_lemu.discovery %>%
  filter(FFGN == 4  & FBGN == 4) 


# Gene Set Enrichment Analysis with ClusterProfiler
#---------------------------------------------------

#we need a ranked list of genes (we can think in using all discovered ones in a similar way than in gene expression we use all not only the differentially expressed ones?)
#determines whether a pre-defined set of genes (ex: those beloging to a specific GO term or KEGG pathway) shows statistically significant, concordant differences between two biological states. In our case the caas number?

# Prepare imput 
# Perform the left join
all4fam_byFGN <- left_join(all4fam_byFGN, gene_lengths, by = "Gene")
cebi_atel_byFGN <- left_join(cebi_atel_byFGN, gene_lengths, by = "Gene")
cebi_lemu_byFGN <- left_join(cebi_lemu_byFGN, gene_lengths, by = "Gene")
lemu_atel_byFGN <- left_join(lemu_atel_byFGN, gene_lengths, by = "Gene")
cerco_atel_byFGN <- left_join(cerco_atel_byFGN, gene_lengths, by = "Gene")
cerco_cebi_byFGN <- left_join(cerco_cebi_byFGN, gene_lengths, by = "Gene")
cerco_lemu_byFGN <- left_join(cerco_lemu_byFGN, gene_lengths, by = "Gene")

# save 7 family contrast in a list
contrasts <- list(all4fam = all4fam_byFGN, cebi_atel = cebi_atel_byFGN, cebi_lemu = cebi_lemu_byFGN, lemu_atel = lemu_atel_byFGN, cerco_atel = cerco_atel_byFGN, cerco_cebi = cerco_cebi_byFGN, cerco_lemu = cerco_lemu_byFGN)

# Loop
for (name in names(contrasts)) {
  
  cat(paste(name), "contrast")
  contrast.rank <- contrasts[[name]] %>%
    group_by(Gene) %>%
    summarise(CAAS_N = n(),
              pval_mean = mean(Pvalue, na.rm = TRUE),
              gene_length_mean = mean(Gene_length_msa, na.rm = TRUE)) %>%
    mutate(rank = ((CAAS_N/gene_length_mean)*(-log10(pval_mean))))
  
  # we order genes by rank statistic created
  gene_list <- contrast.rank$rank 
  # name the vector
  names(gene_list) <- as.character(contrast.rank$Gene)
  # omit any NA values 
  gene_list <- na.omit(gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  cat("GO BP GSEA")
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
  print(dim(gse))
  
  gse_sign <- as.data.frame(gse) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(gse_sign))
  
  # Export dataframe
  write.table(gse_sign, file.path(resultsDir, paste0("/functional/GSEA_results/GOBP_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  cat("GO MF GSEA")
  set.seed(123)
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
                 nPermSimple = 1000 )
  print(dim(mfgse))
  mfgse_sign <- as.data.frame(mfgse) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(mfgse_sign))
  
  write.table(mfgse_sign, file.path(resultsDir, paste0("/functional/GSEA_results/GOMF_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  cat("GO CC GSEA")
  set.seed(123)
  ccgse <- gseGO(geneList = gene_list, 
                 ont = "CC", 
                 keyType = "SYMBOL",
                 minGSSize = 15, 
                 maxGSSize = 500,
                 pvalueCutoff = 1, 
                 verbose = FALSE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "BH",
                 seed = TRUE,
                 scoreType = "pos",
                 nPermSimple = 1000 )  # Whether to print detailed information
  
  print(dim(ccgse))
  ccgse_sign <- as.data.frame(ccgse) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(ccgse_sign))
  
  write.table(ccgse_sign, file.path(resultsDir, paste0("/functional/GSEA_results/GOCC_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # Obtain desired types of IDs from the symbol IDs from significant genes
  ann <- bitr(names(gene_list), fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>%
    filter(!duplicated(SYMBOL))
  # Creating a named vector for mapping
  mapping <- setNames(ann$ENTREZID, ann$SYMBOL)
  
  # Replace names in gene_list with Entrez IDs
  names(gene_list) <- mapping[names(gene_list)]
  
  # Handling any NA values that may have resulted from unmapped symbols
  # Remove NA values (unmapped genes)
  gene_list <- gene_list[!is.na(names(gene_list))]
  
  cat("KEGG GSEA")
  set.seed(123)
  kk_gse <- gseKEGG(geneList     = gene_list,
                    organism     = "hsa",
                    minGSSize    = 15, 
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    nPermSimple = 1000, 
                    pAdjustMethod = "BH",
                    verbose      = FALSE,
                    seed = TRUE,
                    scoreType = "pos")  # Whether to print detailed information
  
  kk_gse <- setReadable(kk_gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  print(dim(kk_gse))
  kkgse_sign <- as.data.frame(kk_gse) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(kkgse_sign))
  
  write.table(kkgse_sign, file.path(resultsDir, paste0("/functional/GSEA_results/KEGG_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  cat("Reactome GSEA")
  set.seed(123)
  fgsea_react <- gsePathway(geneList = gene_list, 
                            organism = "human",
                            minGSSize = 15, 
                            maxGSSize = 500, 
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            eps = 0, 
                            nPermSimple = 1000, 
                            seed = TRUE,
                            scoreType = "pos")  # Whether to print detailed information
  
  fgsea_react <- setReadable(fgsea_react, 'org.Hs.eg.db')
  print(dim(fgsea_react))
  react_sign <- as.data.frame(fgsea_react) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(react_sign))
  
  write.table(react_sign, file.path(resultsDir, paste0("/functional/GSEA_results/react_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  cat("DOSE GSEA")
  set.seed(123)
  do_gsea <- gseDO(geneList = gene_list, 
                   minGSSize = 15, 
                   maxGSSize = 500, 
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   nPermSimple = 1000, 
                   seed = TRUE,
                   scoreType = "pos",
                   verbose = FALSE)  # Whether to print detailed information
  
  do_gsea <- setReadable(do_gsea, 'org.Hs.eg.db')
  print(dim(do_gsea))
  do_sign <- as.data.frame(do_gsea) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(do_sign))
  
  write.table(do_sign, file.path(resultsDir, paste0("/functional/GSEA_results/do_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  cat("DGN GSEA")
  set.seed(123)
  dgn_gsea <- gseDGN(geneList = gene_list,
                     minGSSize = 15, 
                     maxGSSize = 500, 
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     nPermSimple = 1000, 
                     seed = TRUE,
                     scoreType = "pos",
                     verbose = FALSE)  # Whether to print detailed information
  
  dgn_gsea <- setReadable(dgn_gsea, 'org.Hs.eg.db')
  print(dim(dgn_gsea))
  dgn_sign <- as.data.frame(dgn_gsea) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(dgn_sign))
  
  write.table(dgn_sign, file.path(resultsDir, paste0("/functional/GSEA_results/dgn_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  cat("HPO GSEA")
  set.seed(123)
  hpo_gsea <- gseHPO(geneList = gene_list,
                     minGSSize = 15, 
                     maxGSSize = 500, 
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     nPermSimple = 1000, 
                     seed = TRUE,
                     scoreType = "pos",
                     verbose = FALSE)  # Whether to print detailed information
  
  hpo_gsea <- setReadable(hpo_gsea, 'org.Hs.eg.db')
  print(dim(hpo_gsea))
  hpo_sign <- as.data.frame(hpo_gsea) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(hpo_sign))
  
  write.table(dgn_sign, file.path(resultsDir, paste0("/functional/GSEA_results/hpo_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  cat("NCG GSEA")
  set.seed(123)
  ncg_gsea <- gseNCG(geneList = gene_list,
                     minGSSize = 15, 
                     maxGSSize = 500, 
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     nPermSimple = 1000, 
                     seed = TRUE,
                     scoreType = "pos",
                     verbose = FALSE)
  ncg_gsea <- setReadable(ncg_gsea, 'org.Hs.eg.db')
  print(dim(ncg_gsea))
  ncg_sign <- as.data.frame(ncg_gsea) %>%
    filter(qvalue <= 0.05)
  cat("Number significant terms"); print(dim(ncg_sign))
  
  write.table(dgn_sign, file.path(resultsDir, paste0("/functional/GSEA_results/NCG_", name,".txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}



#-----------------------------------------------------------------

# Full analysis with cebi lemu contrast
# Perform the left join
cebi_lemu_byFGN <- left_join(cebi_lemu_byFGN, gene_lengths, by = "Gene")

cebi_lemu.rank <- cebi_lemu_byFGN %>%
  group_by(Gene) %>%
  summarise(CAAS_N = n(),
            pval_mean = mean(Pvalue, na.rm = TRUE),
            gene_length_mean = mean(Gene_length_msa, na.rm = TRUE)) %>%
  mutate(rank = ((CAAS_N/gene_length_mean)*(-log10(pval_mean))))

# we order genes by CAAS number x -log(pval)
gene_list <- cebi_lemu.rank$rank 
# name the vector
names(gene_list) <- as.character(cebi_lemu.rank$Gene)
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


clusterProfiler::ridgeplot(gse, showCategory = 10,core_enrichment = TRUE, label_format = 30)
gseaplot(gse, geneSetID = 1, title = gse$Description[1])

set.seed(123)
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
             nPermSimple = 1000 )  # Whether to print detailed information
            
dim(mfgse)
mfgse_sign <- as.data.frame(mfgse) %>%
  filter(qvalue <= 0.05)
dim(mfgse_sign)

mfsign <- mfgse_sign %>%
  mutate(ShortDescription =  str_wrap(Description, width = 30)) %>%
  select(ID, ShortDescription, Description, setSize, NES, qvalue)

# plot significant results as ER barplot
p1 <- ggplot(mfsign[1:7,], aes(x = reorder(ShortDescription, NES), y = NES, fill = -qvalue)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "#377EB8", high = "#dd6765") +
  labs(x = "", y = "Normalized Enrichment Score", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))
ggsave(plot = p1, filename = file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_ERbarplot.png"), dpi = 300, units = "in",width = 8, height = 8)

# for spoting upregulated or downregulated terms (think more usefuf for expression data)
p <- clusterProfiler::ridgeplot(mfgse, showCategory = 7,core_enrichment = TRUE, label_format = 30)
ggsave(plot = p, filename = file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_ridgeplot.png"), dpi = 300,units = "in",width = 8, height = 8)
p2 <- gseaplot(mfgse, geneSetID = 1, title = mfgse$Description[1])
ggsave(plot = p2, filename = file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_enrichmentplot.png"), dpi = 300,units = "in",width = 8, height = 8)
p3 <- dotplot(mfgse, showCategory = 7, font.size = 9)
ggsave(plot = p3, filename = file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_dotplot.png"), dpi = 300, units = "in",width = 6, height = 6)
# Create Emap plot of 10 most significant terms
ekk <- enrichplot::pairwise_termsim(mfgse)
p4 <- emapplot(ekk, showCategory = 7)
ggsave(plot = p4, filename = file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_emapplot.png"), dpi = 300,units = "in",width = 8, height = 8)

# Cenet plot of 10 most significant terms
p5 <- cnetplot(mfgse, showCategory = 7, cex.params = list(category_label = 0.6), categorySize = "pvalue", foldChange = gene_list) + 
  scale_color_gradient2(name = 'Rank score', low = 'lightgreen', high = 'darkgreen')
ggsave(plot = p5, filename = file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_cnetplot.png"), dpi = 300, units = "in",width = 6, height = 6)

# Prepare figure 8
# Read the images
image1 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_ERbarplot.png"))
image2 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_enrichmentplot.png"))
image3 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_emapplot.png"))

# Annotate each image with letters A, B, and C respectively
image1 <- ggdraw() + draw_image(image1)
image2 <- ggdraw() + draw_image(image2)
image3 <- ggdraw() + draw_image(image3)

# Combine the annotated images
combined_plot <- ((image1 | image2) / image3 ) + 
  plot_annotation(tag_levels = 'A')
# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"functional/GSEA_results/GOMFcebi_lemu_all.png"),  dpi = 300)


set.seed(123)
ccgse <- gseGO(geneList = gene_list, 
               ont = "CC", 
               keyType = "SYMBOL",
               minGSSize = 15, 
               maxGSSize = 500,
               pvalueCutoff = 1, 
               verbose = FALSE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "BH",
               seed = TRUE,
               scoreType = "pos",
               nPermSimple = 1000 )  # Whether to print detailed information

dim(ccgse)
ccgse_sign <- as.data.frame(ccgse) %>%
  filter(qvalue <= 0.05)
dim(ccgse_sign)


sign <- ccgse_sign %>%
  mutate(ShortDescription =  str_wrap(Description, width = 30)) %>%
  select(ID, ShortDescription, Description, setSize, NES, qvalue)

# plot significant results as ER barplot
p1 <- ggplot(sign[1:1,], aes(x = reorder(ShortDescription, NES), y = NES, fill = -qvalue)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "#377EB8", high = "#dd6765") +
  labs(x = "", y = "Normalized Enrichment Score", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))
ggsave(plot = p1, filename = file.path(resultsDir,"functional/GSEA_results/GOCC_cebi_ERbarplot.png"), dpi = 300, units = "in",width = 8, height = 8)

# for spoting upregulated or downregulated terms (think more usefuf for expression data)
p <- clusterProfiler::ridgeplot(ccgse, showCategory = 1,core_enrichment = TRUE, label_format = 30)
ggsave(plot = p, filename = file.path(resultsDir,"functional/GSEA_results/GOCC_cebi_ridgeplot.png"), dpi = 300,units = "in",width = 8, height = 8)
p2 <- gseaplot(ccgse, geneSetID = 1, title = ccgse$Description[1])
ggsave(plot = p2, filename = file.path(resultsDir,"functional/GSEA_results/GOCC_cebi_enrichmentplot.png"), dpi = 300,units = "in",width = 8, height = 8)
p3 <- dotplot(ccgse, showCategory = 1, font.size = 9)
ggsave(plot = p3, filename = file.path(resultsDir,"functional/GSEA_results/GOCC_cebi_dotplot.png"), dpi = 300, units = "in",width = 6, height = 6)
# Create Emap plot of 10 most significant terms
ekk <- enrichplot::pairwise_termsim(ccgse)
p4 <- emapplot(ekk, showCategory = 1)
ggsave(plot = p4, filename = file.path(resultsDir,"functional/GSEA_results/GOCC_cebi_emapplot.png"), dpi = 300,units = "in",width = 8, height = 8)

# Cenet plot of 10 most significant terms
p5 <- cnetplot(ccgse, showCategory = 1, cex.params = list(category_label = 0.6), categorySize = "pvalue", foldChange = gene_list) + 
  scale_color_gradient2(name = 'Rank score', low = 'lightgreen', high = 'darkgreen')
ggsave(plot = p5, filename = file.path(resultsDir,"functional/GSEA_results/GOCC_cebi_cnetplot.png"), dpi = 300, units = "in",width = 6, height = 6)

# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(names(gene_list), fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>%
  filter(!duplicated(SYMBOL))
# Creating a named vector for mapping
mapping <- setNames(ann$ENTREZID, ann$SYMBOL)

# Replace names in gene_list with Entrez IDs
names(gene_list) <- mapping[names(gene_list)]

# Handling any NA values that may have resulted from unmapped symbols
# Remove NA values (unmapped genes)
gene_list <- gene_list[!is.na(names(gene_list))]

set.seed(123)
kk_gse <- gseKEGG(geneList     = gene_list,
                  organism     = "hsa",
                  minGSSize    = 15, 
                  maxGSSize = 500,
                  pvalueCutoff = 1,
                  nPermSimple = 1000, 
                  pAdjustMethod = "BH",
                  verbose      = FALSE,
                  seed = TRUE,
                  scoreType = "pos")  # Whether to print detailed information

kk_gse <- setReadable(kk_gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
dim(kk_gse)
kkgse_sign <- as.data.frame(kk_gse) %>%
  filter(qvalue <= 0.05)
dim(kkgse_sign)


kksign <- kkgse_sign %>%
  mutate(ShortDescription =  str_wrap(Description, width = 30)) %>%
  select(ID, ShortDescription, Description, setSize, NES, qvalue)

# plot significant results as ER barplot
p1 <- ggplot(kksign[1:10,], aes(x = reorder(ShortDescription, NES), y = NES, fill = -qvalue)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "#377EB8", high = "#dd6765") +
  labs(x = "", y = "Normalized Enrichment Score", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))
ggsave(plot = p1, filename = file.path(resultsDir,"functional/GSEA_results/kegg_cebi_ERbarplot.png"), dpi = 300, units = "in",width = 8, height = 8)

# for spoting upregulated or downregulated terms (think more usefuf for expression data)
p <- clusterProfiler::ridgeplot(kk_gse, showCategory = 10,core_enrichment = TRUE, label_format = 30)
ggsave(plot = p, filename = file.path(resultsDir,"functional/GSEA_results/kegg_cebi_ridgeplot.png"), dpi = 300,units = "in",width = 8, height = 8)
p2 <- gseaplot(kk_gse, geneSetID = 1, title = kk_gse$Description[1])
ggsave(plot = p2, filename = file.path(resultsDir,"functional/GSEA_results/kegg_cebi_enrichmentplot.png"), dpi = 300,units = "in",width = 9, height = 8)
p3 <- dotplot(kk_gse, showCategory = 10, font.size = 9)
ggsave(plot = p3, filename = file.path(resultsDir,"functional/GSEA_results/kegg_cebi_dotplot.png"), dpi = 300, units = "in",width = 6, height = 6)
# Create Emap plot of 10 most significant terms
ekk <- enrichplot::pairwise_termsim(kk_gse)
p4 <- emapplot(ekk, showCategory = 14)
ggsave(plot = p4, filename = file.path(resultsDir,"functional/GSEA_results/kegg_cebi_emapplot.png"), dpi = 300,units = "in",width = 8, height = 8)

# Cenet plot of 10 most significant terms
p5 <- cnetplot(kk_gse, showCategory = 5, cex.params = list(category_label = 0.6), categorySize = "pvalue", foldChange = gene_list) + scale_color_gradient2(name = 'Rank score', low = 'lightgreen', high = 'darkgreen')
ggsave(plot = p5, filename = file.path(resultsDir,"functional/GSEA_results/kegg_cebi_cnetplot.png"), dpi = 300, units = "in",width = 6, height = 6)

# # Split the 'core_enrichment' column into multiple rows
# shared_genes <- kkgse_sign %>%
#   separate_rows(core_enrichment, sep = '/') %>%
#   group_by(core_enrichment) %>%
#   filter(n_distinct(ID) > 1) %>%
#   ungroup() %>%
#   pull(unique(core_enrichment))     
# 
# # we order genes by CAAS number x -log(pval)
# symbol_list <- cebi_lemu.rank$rank 
# # name the vector
# names(symbol_list) <- as.character(cebi_lemu.rank$Gene)
# # omit any NA values 
# symbol_list <- na.omit(symbol_list)
# # sort the list in decreasing order 
# symbol_list <- sort(symbol_list, decreasing = TRUE)
# 
# # 'symbol_list' is a named vector of rank stat with gene names as names
# shared_geneList <- symbol_list[names(symbol_list) %in% shared_genes]
# 
# # Use the plotEnrich function to plot only the shared genes
# genekitr::plotEnrich(kkgse_sign, 
#                     plot_type = "geneheat", 
#                     show_gene = shared_genes, 
#                     fold_change = shared_geneList)

# Prepare figure 9
# Read the images
image1 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/kegg_cebi_ERbarplot.png"))
image2 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/kegg_cebi_enrichmentplot.png"))
image3 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/kegg_cebi_emapplot.png"))

# Annotate each image with letters A, B, and C respectively
image1 <- ggdraw() + draw_image(image1)
image2 <- ggdraw() + draw_image(image2)
image3 <- ggdraw() + draw_image(image3)

# Combine the annotated images
combined_plot <- ((image1 | image2) / image3 ) + 
  plot_annotation(tag_levels = 'A')
# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"functional/GSEA_results/keggcebi_lemu_all.png"),  dpi = 400)

set.seed(123)
fgsea_react <- gsePathway(geneList = gene_list, 
                          organism = "human",
                          minGSSize = 15, 
                          maxGSSize = 500, 
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          eps = 0, 
                          nPermSimple = 1000, 
                          seed = TRUE,
                          scoreType = "pos")  # Whether to print detailed information
                         
fgsea_react <- setReadable(fgsea_react, 'org.Hs.eg.db')
dim(fgsea_react)
react_sign <- as.data.frame(fgsea_react) %>%
  filter(qvalue <= 0.05)
dim(react_sign)


reactsign <- react_sign %>%
  mutate(ShortDescription =  str_wrap(Description, width = 30)) %>%
  select(ID, ShortDescription, Description, setSize, NES, qvalue)

# plot significant results as ER barplot
p1 <- ggplot(reactsign[1:6,], aes(x = reorder(ShortDescription, NES), y = NES, fill = -qvalue)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low = "#377EB8", high = "#dd6765") +
  labs(x = "", y = "Normalized Enrichment Score", fill = "FDR") +
  theme_classic() +
  theme(axis.text = element_text(size = 11))
ggsave(plot = p1, filename = file.path(resultsDir,"functional/GSEA_results/reactome_cebi_ERbarplot.png"), dpi = 300, units = "in",width = 8, height = 8)

# for spoting upregulated or downregulated terms (think more usefuf for expression data)
p <- clusterProfiler::ridgeplot(fgsea_react, showCategory = 6, core_enrichment = TRUE, label_format = 30)
ggsave(plot = p, filename = file.path(resultsDir,"functional/GSEA_results/reactome_cebi_ridgeplot.png"), dpi = 300,units = "in",width = 8, height = 8)
p2 <- gseaplot(fgsea_react, geneSetID = 1, title = fgsea_react$Description[1])
ggsave(plot = p2, filename = file.path(resultsDir,"functional/GSEA_results/reactome_cebi_enrichmentplot.png"), dpi = 300,units = "in",width = 8, height = 8)
p3 <- dotplot(fgsea_react, showCategory = 6, font.size = 9)
ggsave(plot = p3, filename = file.path(resultsDir,"functional/GSEA_results/reactome_cebi_dotplot.png"), dpi = 300, units = "in",width = 6, height = 6)
# Create Emap plot of 10 most significant terms
ekk <- enrichplot::pairwise_termsim(fgsea_react)
p4 <- emapplot(ekk, showCategory = 6)
ggsave(plot = p4, filename = file.path(resultsDir,"functional/GSEA_results/reactome_cebi_emapplot.png"), dpi = 300,units = "in",width = 8, height = 8)

# Cenet plot of 10 most significant terms
p5 <- cnetplot(fgsea_react, showCategory = 6, cex.params = list(category_label = 0.6), categorySize = "pvalue", foldChange = gene_list) + 
  scale_color_gradient2(name = 'Rank score', low = 'lightgreen', high = 'darkgreen')
ggsave(plot = p5, filename = file.path(resultsDir,"functional/GSEA_results/reactome_cebi_cnetplot.png"), dpi = 300, units = "in",width = 6, height = 6)

# Prepare figure 10
# Read the images
image1 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/reactome_cebi_ERbarplot.png"))
image2 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/reactome_cebi_enrichmentplot.png"))
image3 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/reactome_cebi_emapplot.png"))

# Annotate each image with letters A, B, and C respectively
image1 <- ggdraw() + draw_image(image1)
image2 <- ggdraw() + draw_image(image2)
image3 <- ggdraw() + draw_image(image3)

# Combine the annotated images
combined_plot <- ((image1 | image2) / image3 ) + 
  plot_annotation(tag_levels = 'A')
# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"functional/GSEA_results/reactome_cebi_lemu_all.png"),  dpi = 400)

set.seed(123)
do_gsea <- gseDO(geneList = gene_list, 
                 minGSSize = 15, 
                 maxGSSize = 500, 
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 nPermSimple = 1000, 
                 seed = TRUE,
                 scoreType = "pos",
                 verbose = FALSE)  # Whether to print detailed information
                
do_gsea <- setReadable(do_gsea, 'org.Hs.eg.db')
dim(do_gsea)
do_sign <- as.data.frame(do_gsea) %>%
  filter(qvalue <= 0.05)
dim(do_sign)


set.seed(123)
dgn_gsea <- gseDGN(geneList = gene_list,
       minGSSize = 15, 
       maxGSSize = 500, 
       pvalueCutoff = 1,
       pAdjustMethod = "BH",
       nPermSimple = 1000, 
       seed = TRUE,
       scoreType = "pos",
       verbose = FALSE)  # Whether to print detailed information

dgn_gsea <- setReadable(dgn_gsea, 'org.Hs.eg.db')
dim(dgn_gsea)
dgn_sign <- as.data.frame(dgn_gsea) %>%
  filter(qvalue <= 0.05)
dim(dgn_sign)


# Prepare figure 11
# Read the images
image1 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/GOMF_cebi_cnetplot.png"))
image2 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/kegg_cebi_cnetplot.png"))
image3 <- magick::image_read(file.path(resultsDir,"functional/GSEA_results/reactome_cebi_cnetplot.png"))

# Annotate each image with letters A, B, and C respectively
image1 <- ggdraw() + draw_image(image1)
image2 <- ggdraw() + draw_image(image2)
image3 <- ggdraw() + draw_image(image3)

# Combine the annotated images
combined_plot <- ((image1 | image2) / image3 ) + 
  plot_annotation(tag_levels = 'A')
# Save to a file
ggsave(combined_plot, filename = file.path(resultsDir,"functional/GSEA_results/cnetplot_cebi_lemu_all.png"),  dpi = 500)

set.seed(123)
hpo_gsea <- gseHPO(geneList = gene_list,
                   minGSSize = 15, 
                   maxGSSize = 500, 
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   nPermSimple = 1000, 
                   seed = TRUE,
                   scoreType = "pos",
                   verbose = FALSE)  # Whether to print detailed information

hpo_gsea <- setReadable(hpo_gsea, 'org.Hs.eg.db')
dim(hpo_gsea)
hpo_sign <- as.data.frame(hpo_gsea) %>%
  filter(qvalue <= 0.05)
dim(hpo_sign)

ncg_gsea <- gseNCG(geneList = gene_list,
       minGSSize = 15, 
       maxGSSize = 500, 
       pvalueCutoff = 1,
       pAdjustMethod = "BH",
       nPermSimple = 1000, 
       seed = TRUE,
       scoreType = "pos",
       verbose = FALSE)
ncg_gsea <- setReadable(ncg_gsea, 'org.Hs.eg.db')
dim(ncg_gsea)
ncg_sign <- as.data.frame(ncg_gsea) %>%
  filter(qvalue <= 0.05)
dim(ncg_sign)

packageVersion("DOSE")
# # MeSH enrichment analysis
# library(meshes)
# library(MeSH.Hsa.eg.db)
# emesh2 <- gseMeSH(geneList, 
#                   MeSHDb = "MeSH.Hsa.eg.db", 
#                   ont = "Disease",
#                   minGSSize = 15, 
#                   maxGSSize = 500, 
#                   pvalueCutoff = 1,
#                   pAdjustMethod = "BH",
#                   nPermSimple = 1000, 
#                   seed = 123,
#                   scoreType = "pos",
#                   verbose = FALSE)
