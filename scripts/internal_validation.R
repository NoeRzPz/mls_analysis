#Phylogenetic ANOVA, often associated with the R package phylolm, is a method used to test the association between a continuous trait and one or more predictors, taking into account the phylogenetic relationships among the observations. It's an adaptation of the classic ANOVA that incorporates a phylogenetic tree to control for non-independence among species due to their shared evolutionary history.

#RRPP (Residual Randomization Permutation Procedure) is a method to create a null distribution of a test statistic that is free of phylogenetic signal. It can be used for hypothesis testing in a phylogenetic context and is a way to perform phylogenetic ANOVA when you need to validate a discovered attribute (like amino acids - AAs in your case) within a phylogenetic tree.

# Load packages
library(tidyverse)
library(phytools)
library(ggplot2)

workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out")

# Import df with CAAS and LQ values
data <- read_tsv(file.path(resultsDir,"/functional/sppclassification_results.tab"), col_names = T)
# Load phylogenetic tree 
primates_tree <- read.tree(file.path(dataDir,"science.abn7829_data_s4.nex.tree"))

# Filter genes with more than 5 cases of each level of the factor type_LQ
data <- data %>%
  group_by(Gene) %>%
  filter( sum(type_LQ == "Long") >= 5 & sum(type_LQ == "Short") >= 5) %>%
  ungroup()  # Ungroup if you need to perform further ungrouped operations

# change the reference level of a factor
#gene_data$type_LQ <- relevel(factor(gene_data$type_LQ), ref = "Short")

results <- list()
graphs <- list()
for (gene in unique(data$Gene)) {
  gene_data <- data %>%
    filter(Gene == gene)
  
  # Find the common names between tree labels and trait df
  common_names <- intersect(primates_tree$tip.label, gene_data$label)
  
  # Prune tree, keeping the species from our trait dataframe, drop tips that are NOT in common_names
  pruned_tree <- drop.tip(primates_tree, setdiff(primates_tree$tip.label, common_names))
  
  # Make sure the tree and the data have the same species order
  gene_data <- gene_data[match(pruned_tree$tip.label, gene_data$label), ]
  
  # Run RRPP
  result <- phytools::phylANOVA(tree = pruned_tree, y = gene_data$LQ, x = gene_data$type_LQ, nsim = 1000)
  
  # Store the results
  results[[gene]] <- result
  
  # Plot of trait values by predictor categories
  LQplot <- ggplot(gene_data, aes(x = type_LQ, y = LQ, fill = type_LQ)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Short" = "#CC313D", "Long" = "#189ab4")) +
  labs(x = "CAAS Predictor", y = "LQ Value", title = paste("Phylogenetic ANOVA for", gene)) +
  theme_classic() +
  theme(legend.position = "none")
  
  # Add the plot to our list, named by the gene for easier reference later
  graphs[[gene]] <- LQplot

}

# Assuming 'result' contains a component with a data frame of residuals or some similar measure
# You might need to adjust this to fit the actual structure of the 'result' object

# 'trait' would be the trait of interest, e.g., presence/absence or count of an amino acid
# 'predictor' could be a categorical variable indicating different groups or conditions
# 'nperm' is the number of permutations for the RRPP

pval <- c()
for (gene in names(results)) {
  cat("\nGene:", gene)
  print(results[[gene]]$Pf)
  pval <- c(pval, results[[gene]]$Pf)
}
# Adjust p-values using Benjamini-Hochberg method
padj <- p.adjust(pval, method = "BH")
padj
# print the plot for a specific gene:
print(graphs[["DNAH8"]])
print(graphs[["SRRM5"]])
results[["SRRM5"]]
# Or loop through the list to print all plots
#for (gene in names(graphs)) {
#  print(graphs[[gene]])
#}

# If you want to save all plots to files
for (gene in names(graphs)) {
  ggsave(filename = paste0("plot_", gene, ".png"), plot = graphs[[gene]], width = 10, height = 8)
}

