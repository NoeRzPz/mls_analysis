#Phylogenetic ANOVA, often associated with the R package phylolm, is a method used to test the association between a continuous trait and one or more predictors, taking into account the phylogenetic relationships among the observations. It's an adaptation of the classic ANOVA that incorporates a phylogenetic tree to control for non-independence among species due to their shared evolutionary history.

#It can be used for hypothesis testing in a phylogenetic context and is a way to perform phylogenetic ANOVA when you need to validate a discovered attribute (like amino acids - AAs in our case) within a phylogenetic tree.

# Load packages
library(tidyverse)
library(phytools)
library(ggplot2)
library(nlme)

workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out")

# Import df with CAAS and LQ values
data <- read_tsv(file.path(resultsDir,"/functional/4fam_sppclassification_results.tab"), col_names = T)
data <- read_tsv(file.path(resultsDir,"/functional/cebi_atel_sppclassification_results.tab"), col_names = T)
data <- read_tsv(file.path(resultsDir,"/functional/cebi_lemu_sppclassification_results.tab"), col_names = T)
data <- read_tsv(file.path(resultsDir,"/functional/lemu_atel_sppclassification_results.tab"), col_names = T)
data <- read_tsv(file.path(resultsDir,"/functional/cerco_cebi_sppclassification_results.tab"), col_names = T)
data <- read_tsv(file.path(resultsDir,"/functional/cerco_atel_sppclassification_results.tab"), col_names = T)
data <- read_tsv(file.path(resultsDir,"/functional/cerco_lemu_sppclassification_results.tab"), col_names = T)

# Load phylogenetic tree 
primates_tree <- read.tree(file.path(dataDir,"science.abn7829_data_s4.nex.tree"))

# Find common names between tree labels and trait dataframe
common_names <- intersect(primates_tree$tip.label, unique(data$label))

data <- data %>%
  filter(label %in% common_names)

# Step 1: Identify the groups (hacer para seleccionar por la mediana!)
valid_groups <- data %>%
  group_by(Gene, Position) %>%
  # Filter genes with more than 5 cases of each level of the factor type_LQ
  filter(sum(type_LQ == "Long") >= 5, sum(type_LQ == "Short") >= 5) %>%
  group_by(Gene, Position, type_LQ) %>%
  summarise(mean_LQ = mean(LQ), .groups = 'drop') %>%
  pivot_wider(names_from = type_LQ, values_from = mean_LQ) %>%
  filter(Long > Short) %>%
  select(Gene, Position) # Obtain genes by position that pass the filter

# Step 2: Filter the original data frame
filtered_data <- data %>%
  semi_join(valid_groups, by = c("Gene", "Position"))

# change the reference level of a factor
#gene_data$type_LQ <- relevel(factor(gene_data$type_LQ), ref = "Short")

# Create unique gene-position pairs
unique_gene_positions <- unique(filtered_data[c("Gene", "Position")])

results <- list()
graphs <- list()

# Iterate over each unique gene-position pair
for (i in 1:nrow(unique_gene_positions)) {
  gene <- unique_gene_positions$Gene[i]
  position <- unique_gene_positions$Position[i]
  
  # Filter data for current gene and position
  gene_position_data <- filtered_data %>%
    filter(Gene == gene, Position == position)
  
  # Find common names between tree labels and trait dataframe
  common_names <- intersect(primates_tree$tip.label, gene_position_data$label)
  
  # Prune tree, keeping species from our trait dataframe, drop tips that are NOT in common_names
  pruned_tree <- drop.tip(primates_tree, setdiff(primates_tree$tip.label, common_names))

  # Ensure the tree and data have the same species order
  gene_position_data <- gene_position_data[match(pruned_tree$tip.label, gene_position_data$label), ]
  
  # Make a name vector with the LQ values from or species
  LQ <- gene_position_data$LQ
  names(LQ) <- gene_position_data$label
  # Set a seed for reproducibility
  set.seed(123)

  result <- phytools::phylANOVA(tree = pruned_tree, y = LQ, x = gene_position_data$type_LQ, nsim = 3000, posthoc = FALSE)
  results[[paste(gene, position, sep = "_")]] <- result
  
  # Plot of trait values by predictor categories
  LQplot <- ggplot(gene_position_data, aes(x = type_LQ, y = LQ, fill = type_LQ)) +
    geom_boxplot() +
    geom_point( alpha = 0.5) +  # Add jitter to avoid overplotting
    scale_fill_manual(values = c("Short" = "#CC313D", "Long" = "#189ab4")) +
    labs(x = "CAAS Predictor", y = "LQ Value", title = paste("Phylogenetic permutation for", gene, "at position", position)) +
    theme_classic() +
    theme(legend.position = "none")
  
  # Add the plot to our list, named by the gene_position for easier reference later
  graphs[[paste(gene, position, sep = "_")]] <- LQplot
}

pval <- c()
for (gene_pos in names(results)) {
  cat("\nGene:", gene_pos)
  print(results[[gene_pos]]$Pf)
  pval <- c(pval, results[[gene_pos]]$Pf)
  # Set the name of the last element in pval to gene_pos
  names(pval)[length(pval)] <- gene_pos
}
significant <- pval[which(pval <= 0.05)]
significant 
# Adjust p-values using Benjamini-Hochberg method (not correct because the different positions in the same gene are not independent, think alternative)
padj <- p.adjust(pval, method = "BH")
padj[which(padj <= 0.05)]

#results <- list()
#graphs <- list()
# Iterate over each unique gene-position pair
for (i in 1:nrow(unique_gene_positions)) {
  gene <- unique_gene_positions$Gene[i]
  position <- unique_gene_positions$Position[i]
  
  # Filter data for current gene and position
  gene_position_data <- filtered_data %>%
    filter(Gene == gene, Position == position)
  
  # Find common names between tree labels and trait dataframe
  common_names <- intersect(primates_tree$tip.label, gene_position_data$label)
  
  # Prune tree, keeping species from our trait dataframe, drop tips that are NOT in common_names
  pruned_tree <- drop.tip(primates_tree, setdiff(primates_tree$tip.label, common_names))
  
  # Ensure the tree and data have the same species order
  gene_position_data <- gene_position_data[match(pruned_tree$tip.label, gene_position_data$label), ]
  
  
  # Fit a PGLS model (without the effect of the group)
  #Performing a permutation test while accounting for phylogenetic relationships
  model <- gls(LQ ~ 1, data = gene_position_data, correlation = corBrownian(phy = pruned_tree))
  model_full <- update(model, .~. + type_LQ)  # Model with the group effect
  
  # Calculate the observed difference
  obs_diff <- with(gene_position_data, mean(LQ[type_LQ == "Long"]) - mean(LQ[type_LQ == "Short"]))
  
  # Permute residuals and calculate differences
  n_perm <- 3000
  perm_diffs <- replicate(n_perm, {
    # Set a seed for reproducibility
    set.seed(123)
    # Shuffle residuals
    shuffled_residuals <- residuals(model)[sample(length(residuals(model)))]
    # Add shuffled residuals to predicted values
    shuffled_LQ <- predict(model_full, gene_position_data) + shuffled_residuals
    # Recalculate group means for the shuffled data
    mean(shuffled_LQ[gene_position_data$type_LQ == "Long"]) - mean(shuffled_LQ[gene_position_data$type_LQ == "Short"])
  })
  
  # Calculate p-value
  p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
  # Store the results
  results[[paste(gene, position, sep = "_")]] <- p_value
  
  # Plot of trait values by predictor categories
  LQplot <- ggplot(gene_position_data, aes(x = type_LQ, y = LQ, fill = type_LQ)) +
    geom_boxplot() +
    geom_point( size = 1.5, alpha = 0.5) +  # Add jitter to avoid overplotting
    scale_fill_manual(values = c("Short" = "#CC313D", "Long" = "#189ab4")) +
    labs(x = "CAAS Predictor", y = "LQ Value", title = paste("Phylogenetic permutation for", gene, "at position", position)) +
    theme_classic() +
    theme(legend.position = "none")
  
  # Add the plot to our list, named by the gene_position for easier reference later
  graphs[[paste(gene, position, sep = "_")]] <- LQplot
}
# Assuming 'result' contains a component with a data frame of residuals or some similar measure
# You might need to adjust this to fit the actual structure of the 'result' object

# 'trait' would be the trait of interest, e.g., presence/absence or count of an amino acid
# 'predictor' could be a categorical variable indicating different groups or conditions
# 'nperm' is the number of permutations for the RRPP

pval <- c()
for (gene_pos in names(results)) {
  cat("\nGene:", gene_pos)
  print(results[[gene_pos]])
  pval <- c(pval, results[[gene_pos]])
  # Set the name of the last element in pval to gene_pos
  names(pval)[length(pval)] <- gene_pos
}
# Implement a phylogenetic t-test-like procedure with empirical p-value calculation
#----------------------------------------------------------------------------------
# Implement a phylogenetic U-Mann-Whitney-like procedure with empirical p-value calculation
#Null hypothesis: two samples come from the same population (i.e. have the same median) or, alternatively, observations in one sample tend to be larger than observations in the other
#-------------------------------------------------------------------------------
phylwilcoxtest <- function(tree, x, y, nsim=1000, data) {
  # Ensure that tree is a phylo object
  if (!inherits(tree, "phylo")) {
    stop("tree should be an object of class \"phylo\".")
  }
  
  # Ensure that x and y are column names in data and not separate vectors
  if (!(x %in% names(data)) || !(y %in% names(data))) {
    stop("x and y must be column names in data.")
  }
  
  # Extract x and y values from data, ordered by tree$tip.label
  x_values <- data[[x]]
  y_values <- data[[y]]
  
  # Convert x to factor
  x_values <- as.factor(x_values)
  # Check if factor has exactly 2 levels
  if (length(levels(x_values)) != 2) {
    stop("Factor does not have 2 levels.")
  }
  
  # Compute BM rate for y
  sig2 <- mean(pic(y_values, multi2di(tree, random = FALSE))^2)
  
  # Calculate the observed difference
  result <- wilcox.test(y_values ~ x_values, alternative = "greater", paired = FALSE, conf.int = 0.95)
  
  W.obs <- result$statistic
  
  # Simulate
  sims <- fastBM(tree, sig2 = sig2, nsim = (nsim - 1))
  W.null <- vector()
  W.null[1] <- W.obs
  
  for (i in 2:nsim) {
    # Recalculate group medians for the simulated data
    simresult <- wilcox.test(sims[,i-1] ~ x_values, alternative = "greater", paired = FALSE, conf.int = 0.95)
    W.sim <- simresult$statistic
    
    W.null[i] <- W.sim
  }
  
  # Calculate p-value
  P.W <- sum(W.null >= W.obs) / nsim
  
  # Return result
  obj <- list(Wobs = W.obs, Pval = P.W)
  class(obj) <- "phylt-test"
  obj
}

results <- list()
graphs <- list()
# Iterate over each unique gene-position pair
for (i in 1:nrow(unique_gene_positions)) {
  gene <- unique_gene_positions$Gene[i]
  position <- unique_gene_positions$Position[i]
  
  # Filter data for current gene and position
  gene_position_data <- filtered_data %>%
    filter(Gene == gene, Position == position)
  
  # Find common names between tree labels and trait dataframe
  common_names <- intersect(primates_tree$tip.label, gene_position_data$label)
  
  # Prune tree, keeping species from our trait dataframe, drop tips that are NOT in common_names
  pruned_tree <- drop.tip(primates_tree, setdiff(primates_tree$tip.label, common_names))
  
  # Ensure the tree and data have the same species order
  gene_position_data <- gene_position_data[match(pruned_tree$tip.label, gene_position_data$label), ]
  # Make a name vector with the LQ values from or species
  LQ <- gene_position_data$LQ
  names(LQ) <- gene_position_data$label
  # Set a seed for reproducibility
  set.seed(123)
  
  result <- phylttest(tree = pruned_tree, y = "LQ", x = "type_LQ", nsim = 3000, data = gene_position_data)
  
  # Store the results
  results[[paste(gene, position, sep = "_")]] <- result
  
  print(paste(paste(gene, position, sep = "_"), "processed"))
  
  # Plot of trait values by predictor categories
  LQplot <- ggplot(gene_position_data, aes(x = type_LQ, y = LQ, fill = type_LQ)) +
    geom_boxplot() +
    geom_point(size = 1.5, alpha = 0.5) +  # Add jitter to avoid overplotting
    scale_fill_manual(values = c("Short" = "#CC313D", "Long" = "#189ab4")) +
    labs(x = "CAAS Predictor", y = "LQ Value", title = paste("Phylogenetic permutation for", gene, "at position", position)) +
    theme_classic() +
    theme(legend.position = "none")
  
  # Add the plot to our list, named by the gene_position for easier reference later
  graphs[[paste(gene, position, sep = "_")]] <- LQplot
}

pval <- c()
for (gene_pos in names(results)) {
  #cat("\nGene:", gene_pos)
  #print(results[[gene_pos]]$Pval)
  pval <- c(pval, results[[gene_pos]]$Pval)
  # Set the name of the last element in pval to gene_pos
  names(pval)[length(pval)] <- gene_pos
  write.table(cbind(names(pval), pval), file.path(resultsDir,"/functional/Internal_validation/phylwilcoxon_pvalues.txt"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}
significant <- pval[which(pval <= 0.05)]
significant 
# Adjust p-values using Benjamini-Hochberg method (not correct because the different positions in the same gene are not independent, think alternative)
padj <- p.adjust(pval, method = "BH")
padj[which(padj <= 0.05)]
# Adjust p-values for FDR (these are the q-values)
q_values <- p.adjust(pval, method = "fdr")
# Display q-values
q_values[which(q_values <= 0.05)]

# loop through the list to print significant plots
for (gene_pos in names(significant)) {
  print(graphs[[gene_pos]])
  # save all significant plots to files
  ggsave(filename = file.path(resultsDir,paste0("/functional/Internal_validation/boxplot_", gene_pos, ".png")), plot = graphs[[gene_pos]], width = 10, height = 8)
}
hist(pval, breaks = 50)
