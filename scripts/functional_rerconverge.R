# Functional enrichment detection methods to find functionally-related genomic elements experiencing convergent evolutionary rates as a group

# Load package
library(RERconverge)
library(ggplot2)
library(ggrepel)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Import inputs: 
# output from RERconverge correlation functions 
res <- read.csv(file.path(resultsDir,"rerconverge/continuostrait_cor_res.csv"), row.names = 1) 
rer_genes <- read.csv(file.path(resultsDir,"rerconverge/signif_lq_genes_005.csv"), row.names = 1)

# Import Pathway Annotations
c2annots <- read.gmt("data/c2.all.v2023.2.Hs.symbols.gmt")
hannots <- read.gmt("data/h.all.v2023.2.Hs.symbols.gmt")
c4annots <- read.gmt("data/c4.all.v2023.2.Hs.symbols.gmt")
c5annots <- read.gmt("data/c5.all.v2023.2.Hs.symbols.gmt")
c6annots <- read.gmt("data/c6.all.v2023.2.Hs.symbols.gmt")
c7annots <- read.gmt("data/c7.all.v2023.2.Hs.symbols.gmt")
c8annots <- read.gmt("data/c8.all.v2023.2.Hs.symbols.gmt")

# Create a ranked gene list
stats <- getStat(res)

# RERconverge enrichment functions expect pathways to be in named pathway-group lists contained within a list 
annotations <- list(c2 = c2annots, h = hannots, c4 = c4annots, c5 = c5annots, c6 = c6annots, c7 = c7annots, c8 = c8annots)

# Define empty results list
enrichmen_res <- list()

for (name in names(annotations)) {
  ann <- paste0(name,"annotlist")
  ann <- list(annotations[[name]])
  names(ann) <- "MSigDBpathways"
  
  # Calculate Enrichment
  # Use Wilcoxon Rank-Sum Test on a list of genes ranked based on their correlation statistics
  enrichment <- fastwilcoxGMTall(stats, ann, outputGeneVals = T, num.g = 10)
  
  # Add result to list
  enrichmen_res[[name]] <- enrichment
}


# Sort the data by AdjustedPValue and select significant results
c2signif_res <- enrichmen_res[["c2"]]$MSigDBpathways %>%
  arrange(p.adj) %>%
  filter(p.adj <= 0.05) 
# No significant enriched terms

hsignif_res <- enrichmen_res[["h"]]$MSigDBpathways %>%
  arrange(p.adj) %>%
  filter(p.adj <= 0.05) 
# No significant enriched terms

#Computational gene sets defined by mining large collections of cancer-oriented expression data.
c4signif_res <- enrichmen_res[["c4"]]$MSigDBpathways %>%
  arrange(p.adj) %>%
  filter(p.adj <= 0.05) 

# Create a bar plot
ggplot(c4signif_res, aes(x = reorder(row.names(c4signif_res), p.adj), y = stat, fill = p.adj)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene Set", y = "Enrichment Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient(low = "blue", high = "red") +
  ggtitle("Enrichment Analysis Results") +
  scale_y_continuous(expand = c(0, 0))
ggsave("out/rerconverge/barplot_c4.png", dpi = 300)

# Create a scatterplot
top_terms <- head(c4signif_res, 5)
ggplot(enrichmen_res[["c4"]]$MSigDBpathways, aes(x = stat, y =  -log10(p.adj), color = p.adj <= 0.05)) +
  geom_point(aes(size = num.genes), alpha = 0.5) +
  geom_text_repel(data = top_terms, aes(label = row.names(top_terms)), 
                  box.padding = 0.5, point.padding = 0.2, 
                  max.overlaps = Inf, nudge_x = 0.1, nudge_y = 0.1,
                  size = 3)  +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(x = "Statistic (stat)", y = "-log10(Adjusted P-Value)") +
  scale_size_continuous(range = c(2, 10)) +  # Adjust the size range as needed
  theme_minimal() +
  guides(color = "none") +
  ggtitle("Enrichment Analysis Results")
ggsave("out/rerconverge/scatterplot_c4.png", dpi = 300)

#Ontology gene sets: GO and HPO
c5signif_res <- enrichmen_res[["c5"]]$MSigDBpathways %>%
  arrange(p.adj) %>%
  filter(p.adj <= 0.05) 
# No significant enriched terms

# Oncogenic signature gene sets,signatures of cellular pathways which are often dis-regulated in cancer
c6signif_res <- enrichmen_res[["c6"]]$MSigDBpathways %>%
  arrange(p.adj) %>%
  filter(p.adj <= 0.05) 

#Immunology signature gene sets: Gene sets that represent cell states and perturbations within the immune system
c7signif_res <- enrichmen_res[["c7"]]$MSigDBpathways %>%
  arrange(p.adj) %>%
  filter(p.adj <= 0.05) 

# Visualize data 
# Create a bar plot
ggplot(c7signif_res, aes(x = reorder(row.names(c7signif_res), p.adj), y = stat, fill = p.adj)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene Set", y = "Enrichment Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient(low = "blue", high = "red") +
  ggtitle("Enrichment Analysis Results") +
  scale_y_continuous(expand = c(0, 0))
ggsave("out/rerconverge/scatterplot_c7.png", dpi = 300)

# Create a scatterplot
top_terms <- head(c7signif_res, 5)
ggplot(enrichmen_res[["c7"]]$MSigDBpathways, aes(x = stat, y =  -log10(p.adj), color = p.adj <= 0.05)) +
  geom_point(aes(size = num.genes), alpha = 0.5) +
  geom_text_repel(data = top_terms, aes(label = row.names(top_terms)), 
                  box.padding = 0.5, point.padding = 0.2, 
                  max.overlaps = Inf, nudge_x = 0.1, nudge_y = 0.1,
                  size = 3)  +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  labs(x = "Statistic (stat)", y = "-log10(Adjusted P-Value)") +
  scale_size_continuous(range = c(2, 10)) +  # Adjust the size range as needed
  theme_minimal() +
  guides(color = "none") +
  ggtitle("Enrichment Analysis Results")
ggsave("out/rerconverge/scatterplot_c7.png", dpi = 300)

#cell type signature gene sets,contain curated cluster markers for cell types identified in single-cell sequencing studies of human tissue
c8signif_res <- enrichmen_res[["c8"]]$MSigDBpathways %>%
  arrange(p.adj) %>%
  filter(p.adj <= 0.05) 


# ORA analyses
#-------------

# Obtain desired types of IDs from the symbol IDs from significant genes
ann <- bitr(unique(rer_genes), fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

# Do the same for all genes studied
background_ann <- bitr(unique(rownames(res)), fromType = "SYMBOL", toType =  c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) %>% 
  filter(!duplicated(SYMBOL))

ego <- enrichGO(gene          = unique(ann$ENTREZID),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "BP",
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff  = 1, 
                qvalueCutoff  = 0.05, 
                readable = TRUE, 
                universe = unique(background_ann$ENTREZID))

dim(ego)
# No enriched GO BP terms for nominal significant RERconverge genes, with minGSSize =15 the results were the same in all cases



ego <- enrichGO(gene          = unique(ann$ENTREZID),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "MF",
                pAdjustMethod = "BH",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff  = 1, 
                qvalueCutoff  = 0.05, 
                readable = TRUE, 
                universe = unique(background_ann$ENTREZID))

dim(ego)
# No enriched GO MF terms for nominal significant RERconverge genes

set.seed(123)
kk <- enrichKEGG(gene         = unique(ann$ENTREZID),
                 organism     = 'hsa',
                 universe = as.character(unique(background_ann$ENTREZID)),
                 pvalueCutoff = 1,
                 qvalueCutoff  = 0.05, 
                 minGSSize     = 10,
                 maxGSSize     = 500)
dim(kk)
# No enriched KEGG terms for nominal significant RERconverge genes

kk <- enrichMKEGG(gene = unique(ann$ENTREZID),
                  universe = unique(background_ann$ENTREZID),
                  organism = 'hsa',
                  pvalueCutoff = 1,
                  qvalueCutoff = 0.05,
                  minGSSize     = 10,
                  maxGSSize     = 500)
dim(kk)
# No enriched KEGG modules terms for nominal significant RERconverge genes

kk <- enrichPathway(gene = unique(ann$ENTREZID),
                    universe = unique(background_ann$ENTREZID), 
                    pvalueCutoff = 1,
                    qvalueCutoff = 0.05,
                    minGSSize     = 10,
                    maxGSSize     = 500,
                    readable = TRUE)
dim(kk)
# No enriched REACTOME terms for nominal significant RERconverge genes


kk <- enrichDO(gene          = unique(ann$ENTREZID),
               ont           = "DO",
               pvalueCutoff  = 1,
               pAdjustMethod = "BH",
               universe      = unique(background_ann$ENTREZID),
               minGSSize     = 10,
               maxGSSize     = 500,
               qvalueCutoff  = 0.05,
               readable      = TRUE)
dim(kk)
# No enriched DOSE terms for nominal significant RERconverge genes



kk <- enrichNCG(unique(ann$ENTREZID),
                pvalueCutoff  = 1,
                pAdjustMethod = "BH",
                universe      = unique(background_ann$ENTREZID),
                minGSSize     = 10,
                maxGSSize     = 500,
                qvalueCutoff  = 0.05,
                readable = TRUE) 
dim(kk)
# No enriched DOSE terms for nominal significant RERconverge genes

kk <- enrichDGN(unique(ann$ENTREZID),
                pAdjustMethod = "BH",
                universe      = unique(background_ann$ENTREZID),
                minGSSize     = 10,
                maxGSSize     = 500,
                pvalueCutoff  = 1,
                qvalueCutoff  = 0.05,
                readable = TRUE) 
dim(kk)
# No enriched DisGeNET terms for nominal significant RERconverge genes

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C4") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C8") %>% 
  dplyr::select(gs_name, entrez_gene)

kk <- enricher(unique(ann$ENTREZID),
               pAdjustMethod = "BH", 
               universe      = unique(background_ann$ENTREZID),
               minGSSize     = 10,
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
               minGSSize     = 10,
               maxGSSize     = 500,
               pvalueCutoff  = 1,
               qvalueCutoff  = 0.05,
               TERM2GENE = m_t2g)
dim(kk)
# No enriched msigdb terms for nominal significant RERconverge genes


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
# No enriched ageing terms for nominal significant RERconverge genes

