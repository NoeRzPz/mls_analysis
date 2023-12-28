# Load library
library(VennDiagram)
library(RColorBrewer)
library(readr)
library(tidyverse)
library(ggplot2)
library(ggvenn)
library(UpSetR)
library(ComplexUpset)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load input data
# ferre validated by RRPP 
farre.discovery <- read_delim(file.path(dataDir,"/ferre_validatedRRPP.txt"),col_names = FALSE ,
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)$X1

farre.discovery <- as.vector(unique(farre.discovery))

# ferre validated by RRPP 
gerard.ILS <- read_delim(file.path(dataDir,"/gerard_increaseLS.txt"),col_names = FALSE ,
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)$X1

gerard.ILS <- as.vector(unique(gerard.ILS))
 
# Load 4 family contrast discovery data frame 
fam4caas.filtered <- read_delim(file.path(resultsDir,"/functional/genes_byFGN.tab"), col_names = FALSE ,
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)$X1
fam4caas.filtered <- as.vector(unique(fam4caas.filtered))

# Load CAAS input data from 2 family contrasts 
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



rer_genes <- read.csv(file.path(resultsDir,"rerconverge/signif_lq_genes_005.csv"), row.names = 1)
genes_shared_all2fam <- read_tsv(file.path(resultsDir,"functional/shared_genesbyall.tsv"), col_names = F)$X1
cebi_lemu.bootstrap <- read_tsv(file.path(resultsDir,"genes_cebi_lemut.tsv"), col_names = F)$X1

# Obtain gene lists without duplicates
genes_cebi_lemu <- unique(cebi_lemu$Gene)

rer_genes <- unique(rownames(rer_genes))

genes_shared_all2fam <- unique(genes_shared_all2fam)
# Obtain desired types of IDs from the symbol IDs from common genes to all contrasts
#ann <- bitr(genes_shared_all2fam, fromType = "SYMBOL", toType =  c("SYMBOL","ENSEMBL","ENTREZID"), OrgDb = org.Hs.eg.db) %>% 
#  filter(!duplicated(SYMBOL))

# Export IDs genes with significant bootstrap aa substitutions for using it in S-LDSC
#write.table(as.vector(ann$ENSEMBL), file = file.path(dataDir,"common_genes.tsv"), sep = "\t",col.names = FALSE, row.names = FALSE, quote = FALSE)

# discovery genes enrichment
lists <- list(`Discovered genes` = genes_cebi_lemu, `Validated genes in mammals` = farre.discovery, `Genes from Gerard et al.` = gerard.ILS )
ggvenn(lists, stroke_size = 0.5, set_name_size = 4,show_percentage = FALSE)

lists2 <- list(`Discovered genes` = genes_cebi_lemu, `Validated genes in mammals` = farre.discovery)
v <- ggvenn(lists2, stroke_size = 0.5, set_name_size = 6, show_percentage = TRUE,text_size = 6) 

# Add annotations
v <- v + annotate("text", x = 0, y = 0, label = "CUBN", size = 5, color = "black", vjust = 4) +
  annotate("text", x = 0, y = 0, label = "FREM2", size = 5, color = "black", vjust = 5.5)

ggsave(v,filename = file.path(resultsDir,paste0("/functional/external_validation/venn_lemu_cebi_mammal.png")),  dpi =  300)

lists3 <- list(`Discovered genes` = genes_cebi_lemu,  `Genes from Gerard et al.` = gerard.ILS )
v <- ggvenn(lists3, stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE)
ggsave(v,filename = file.path(resultsDir,paste0("/functional/external_validation/venn_lemu_cebi_gerard.png")),  dpi =  300)

# Hypergeometric test
# Total number of genes in the genome (or your universe of genes)
total_genes <- 16133
results <- list()
for (contrast in names(lists)[c(1,3)]) {
  # Calculate the overlap
  overlap <- length(intersect(lists[[contrast]], farre.discovery))
  
  # Create a contingency table
  contingency_table <- matrix(c(
    overlap, length(lists[[contrast]]) - overlap,  # Overlap and unique in df1
    length(farre.discovery) - overlap, total_genes - length(lists[[contrast]]) - length(farre.discovery) + overlap  # Unique in df2 and neither
  ), nrow = 2)
  
  # Perform Fisher's Exact Test
  # test for over-representation in the overlap,to find out if there's a higher-than-expected overlap between two sets of genes,
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  
  results[[contrast]] <- fisher_result
}
results

total_genes <- 16133
results <- list()
for (contrast in names(lists)[c(1,2)]) {
  # Calculate the overlap
  overlap <- length(intersect(lists[[contrast]], gerard.ILS))
  
  # Create a contingency table
  contingency_table <- matrix(c(
    overlap, length(lists[[contrast]]) - overlap,  # Overlap and unique in df1
    length(gerard.ILS) - overlap, total_genes - length(lists[[contrast]]) - length(gerard.ILS) + overlap  # Unique in df2 and neither
  ), nrow = 2)
  
  # Perform Fisher's Exact Test
  # test for over-representation in the overlap,to find out if there's a higher-than-expected overlap between two sets of genes,
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  
  results[[contrast]] <- fisher_result
}
results
# Find common elements between 2 sets
common_genes <- intersect(genes_cebi_lemu, farre.discovery)
common_genes
#how many of the overlaping are internally validated
intersect(common_genes, cebi_lemu_validated)
# shared by 3 datasets
shared_elements <- Reduce(intersect, lists)
shared_elements

# Find common elements between 2 primates sets
common_genes <- intersect(genes_cebi_lemu, gerard.ILS)
common_genes

# bootstrap genes enrichment
lists <- list(`Primates genes` = unique(cebi_lemu.bootstrap), `Mammals genes Farre` = farre.discovery, `Primates Gerard` = gerard.ILS )
ggvenn(lists, stroke_size = 0.5, set_name_size = 4,show_percentage = FALSE)

# Hypergeometric test
# Total number of genes in the genome (or your universe of genes)
total_genes <- 16133
results <- list()
for (contrast in names(lists)[c(1,3)]) {
  # Calculate the overlap
  overlap <- length(intersect(lists[[contrast]], farre.discovery))
  
  # Create a contingency table
  contingency_table <- matrix(c(
    overlap, length(lists[[contrast]]) - overlap,  # Overlap and unique in df1
    length(farre.discovery) - overlap, total_genes - length(lists[[contrast]]) - length(farre.discovery) + overlap  # Unique in df2 and neither
  ), nrow = 2)
  
  # Perform Fisher's Exact Test
  # test for over-representation in the overlap,to find out if there's a higher-than-expected overlap between two sets of genes,
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  
  results[[contrast]] <- fisher_result
}
results

# Assume you have two sets of genes (for example, differentially expressed genes from two different dataframes)
pairs <- list("CAAStools significant genes" = genes_cebi_atel, "RERconverge significant genes" = rer_genes)
pairs <- list("CAAStools significant genes" = genes_lemu_atel, "RERconverge significant genes" = rer_genes)

pairs <- list("CAAStools significant genes" = genes_cebi_lemu, "RERconverge significant genes" = rer_genes)
ggvenn(pairs, stroke_size = 0.5, set_name_size = 4,show_percentage = FALSE)
ggsave("out/rerconverge/venn_cebi.lemu.png", dpi = 300)

pair <- list("CAAStools significant bootstrap genes" = cebi_lemu.bootstrap, "RERconverge significant genes" = rer_genes)
ggvenn(pair, stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE)
ggsave("out/rerconverge/venn_cebi.lemubootstrap.png", dpi = 300)

pairs <- list("CAAStools significant genes in common" = genes_shared_all2fam, "RERconverge significant genes" = rer_genes)
ggvenn(pairs, stroke_size = 0.5, set_name_size = 4, fill_color = c("green","yellow"),show_percentage = FALSE)
ggsave("out/rerconverge/venn_common_genes.png", dpi = 300)

# Total number of genes in the genome (or your universe of genes)
total_genes <- 16133

results <- list()
for (contrast in names(lists)) {
  # Calculate the overlap
  overlap <- length(intersect(lists[[contrast]], rer_genes))
  
  # Create a contingency table
  contingency_table <- matrix(c(
    overlap, length(lists[[contrast]]) - overlap,  # Overlap and unique in df1
    length(rer_genes) - overlap, total_genes - length(lists[[contrast]]) - length(rer_genes) + overlap  # Unique in df2 and neither
  ), nrow = 2)
  
  # Perform Fisher's Exact Test
  # test for over-representation in the overlap,to find out if there's a higher-than-expected overlap between two sets of genes,
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  
  results[[contrast]] <- fisher_result
}
# Enrichment for bootstrap genes
for (contrast in names(pairs)) {
  # Calculate the overlap
  overlap <- length(intersect(pairs[[contrast]], rer_genes))
  
  # Create a contingency table
  contingency_table <- matrix(c(
    overlap, length(pairs[[contrast]]) - overlap,  # Overlap and unique in df1
    length(rer_genes) - overlap, total_genes - length(pairs[[contrast]]) - length(rer_genes) + overlap  # Unique in df2 and neither
  ), nrow = 2)
  
  # Perform Fisher's Exact Test
  # test for over-representation in the overlap,to find out if there's a higher-than-expected overlap between two sets of genes,
  fisher_result <- fisher.test(contingency_table, alternative = "greater")
  
  results[[contrast]] <- fisher_result
}
results
# Print the p-value
#fisher_result$p.value
results <- list()
for (contrast in names(lists)) {
  # Calculate the overlap
  overlap <- length(intersect(lists[[contrast]], rer_genes))
  
  #  using a hypergeometric test
  M <- 16133 # total number of genes tested
  N <- length(lists[[contrast]]) # number of genes in the first set
  n <- length(rer_genes) # number of genes in the second set
  k <- length(overlap) # number of common genes
  
  p_value <- phyper(k - 1, N, M - N, n, lower.tail = FALSE)
  results[[contrast]] <- p_value
}
results

# Cebi_lemu internaly validaded genes
cebi_lemu_validated <- c("CUBN", "C6orf89", "CA6", "CARS"," FADS2", "HEMGN", "HIVEP1", "KIAA1671","KRT84", "MFHAS1", "RSPH3", "SETD2", "TMX4",  "ARID3C", "ARMT1", "ASGR1", "AURKA", "C1orf87", "CCDC85A", "CTRC", "EPN3", "FREM2","FUBP3", "GBGT1", "GLI4", "GPR107", "IGF2R", "MTF1", "RHBG", "TACR2", "TCN1", "TG", "TTLL7")

lists <- list(Primates_CAAS = unique(cebi_lemu_validated), Mammals_CAAS_Ferre = farre.discovery, Primates_Gerard = gerard.ILS )

ggvenn(lists, stroke_size = 0.5, set_name_size = 4,show_percentage = FALSE)


# Load internally validated CAAS
all4fam <- read.table(file.path(resultsDir,"/functional/Internal_validation/4fam/phylwilcoxon_pvalues.txt"), row.names = 1)

cebi_atel <- read.table(file.path(resultsDir,"/functional/Internal_validation/cebi_atel/phylwilcoxon_pvalues.txt"), row.names = 1)

cebi_lemu <- read.table(file.path(resultsDir,"/functional/Internal_validation/cebi_lemu/phylwilcoxon_pvalues.txt"), row.names = 1)

lemu_atel <- read.table(file.path(resultsDir,"/functional/Internal_validation/lemu_atel/phylwilcoxon_pvalues.txt"), row.names = 1)

cerco_atel <- read.table(file.path(resultsDir,"/functional/Internal_validation/cerco_atel/phylwilcoxon_pvalues.txt"), row.names = 1)

cerco_cebi <- read.table(file.path(resultsDir,"/functional/Internal_validation/cerco_cebi/phylwilcoxon_pvalues.txt"), row.names = 1)

cerco_lemu <- read.table(file.path(resultsDir,"/functional/Internal_validation/cerco_lemu/phylwilcoxon_pvalues.txt"), row.names = 1)

# Check duplicates
caas_all4fam <- all4fam %>%
  filter(pval <= 0.05) %>%
  tibble::rownames_to_column(var = "caas_name") %>%
  pull(caas_name)
  
caas_cebi_atel <- cebi_atel %>%
  filter(pval <= 0.05) %>%
  tibble::rownames_to_column(var = "caas_name") %>%
  pull(caas_name)

caas_cebi_lemu <- cebi_lemu %>%
  filter(pval <= 0.05) %>%
  tibble::rownames_to_column(var = "caas_name") %>%
  pull(caas_name)

caas_lemu_atel <- lemu_atel %>%
  filter(pval <= 0.05) %>%
  tibble::rownames_to_column(var = "caas_name") %>%
  pull(caas_name)

caas_cerco_atel <- cerco_atel %>%
  filter(pval <= 0.05) %>%
  tibble::rownames_to_column(var = "caas_name") %>%
  pull(caas_name)

caas_cerco_cebi <- cerco_cebi %>%
  filter(pval <= 0.05) %>%
  tibble::rownames_to_column(var = "caas_name") %>%
  pull(caas_name)

caas_cerco_lemu <- cerco_lemu %>%
  filter(pval <= 0.05) %>%
  tibble::rownames_to_column(var = "caas_name") %>%
  pull(caas_name)

# Create a list of your lists
lists <- list(all4fam = caas_all4fam, Cebidae.Atelidae = caas_cebi_atel, Cebidae.Lemuridae = caas_cebi_lemu, Lemuridae.Atelidae = caas_lemu_atel, Cercopithecidae.Atelidae = caas_cerco_atel, Cercopithecidae.Cebidae = caas_cerco_cebi, Cercopithecidae.Lemuridae = caas_cerco_lemu)

# Convert to binary membership data
upset_data <- fromList(lists)

# Create a named vector of colors for the sets with RColorBrewer 
set_colors <- brewer.pal(7, "Set3")
names(set_colors) <- names(lists)

# Generate the UpSet plot with customized colors
upset(upset_data, sets = names(lists), sets.bar.color = set_colors)

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
elements_in_3_sets <- element_counts %>%
  filter(count == 3) %>%
  pull(element)

# Print the elements
print(elements_in_3_sets)