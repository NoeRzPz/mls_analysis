#Phylogenetic ANOVA, often associated with the R package phylolm, is a method used to test the association between a continuous trait and one or more predictors, taking into account the phylogenetic relationships among the observations. It's an adaptation of the classic ANOVA that incorporates a phylogenetic tree to control for non-independence among species due to their shared evolutionary history.

#It can be used for hypothesis testing in a phylogenetic context and is a way to perform phylogenetic ANOVA when you need to validate a discovered attribute (like amino acids - AAs in our case) within a phylogenetic tree.

# Load packages
library(readr)
library(tidyverse)
library(phytools)
library(ggplot2)
library(nlme)
library(reticulate)
library(processx)
library(plotly)
library(cowplot)
library(patchwork)
library(magick)
library(rentrez)

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
  summarise(mean_LQ = median(LQ, na.rm = TRUE), .groups = 'drop') %>% # mean
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

# pval <- c()
# for (gene_pos in names(results)) {
#   cat("\nGene:", gene_pos)
#   print(results[[gene_pos]])
#   pval <- c(pval, results[[gene_pos]])
#   # Set the name of the last element in pval to gene_pos
#   names(pval)[length(pval)] <- gene_pos
# }
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
  
  result <- phylwilcoxtest(tree = pruned_tree, y = "LQ", x = "type_LQ", nsim = 3000, data = gene_position_data)
  
  # Store the results
  results[[paste(gene, position, sep = "_")]] <- result
  
  print(paste(paste(gene, position, sep = "_"), "processed"))
  
  # Plot of trait values by predictor categories
  LQplot <- ggplot(gene_position_data, aes(x = type_LQ, y = LQ, fill = type_LQ)) +
    geom_boxplot() +
    geom_point(size = 1.5, alpha = 0.5) +  # Add jitter to avoid overplotting
    scale_fill_manual(values = c("Short" = "#CC313D", "Long" = "#189ab4")) +
    labs(x = "CAAS Predictor", y = "LQ Value", title = paste( gene, " gene at position", position)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 14),  
          axis.title.y = element_text(size = 14),  
          plot.title = element_text(size = 16))
  
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

#---------------
# Create barplots to spot chromosomes where internally validated genes fall

# List of internaly validated genes you want to query
all4fam_validated <- read_delim(file.path(resultsDir,"/functional/Internal_validation/all4fam/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cebi_atel_validated <- read_delim(file.path(resultsDir,"/functional/Internal_validation/cebi_atel/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cebi_lemu_validated <- read_delim(file.path(resultsDir,"/functional/Internal_validation/cebi_lemu/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cerco_atel_validated <- read_delim(file.path(resultsDir,"/functional/Internal_validation/cerco_atel/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cerco_cebi_validated <- read_delim(file.path(resultsDir,"/functional/Internal_validation/cerco_cebi/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cerco_lemu_validated <- read_delim(file.path(resultsDir,"/functional/Internal_validation/cerco_lemu/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
lemu_atel_validated <- read_delim(file.path(resultsDir,"/functional/Internal_validation/lemu_atel/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)

# Define function to fetch chromosome numbers
get_chromosome_number <- function(gene_symbol) {
  gene_search <- entrez_search(db = "gene", term = paste(gene_symbol, "[Gene Name] AND Homo sapiens[Organism]"))
  
  if (length(gene_search$ids) > 0) {
    gene_summary <- entrez_summary(db = "gene", id = gene_search$ids[1])
    return(gene_summary$chromosome)
  } else {
    return(NA)  # Return NA if gene is not found
  }
}
# save 7 family contrast in a list
contrasts <- list(all4fam = all4fam_validated, cebi_atel = cebi_atel_validated, cebi_lemu = cebi_lemu_validated, cerco_atel = cerco_atel_validated, cerco_cebi = cerco_cebi_validated, cerco_lemu = cerco_lemu_validated, lemu_atel = lemu_atel_validated)

# Loop
for (name in names(contrasts)) {
  cat(paste(name), "contrast")
  genes <- contrasts[[name]]  %>%
    filter(pval <= 0.05) %>% # Filter rows where pval is less than or equal to 0.05
    mutate(gene = sapply(strsplit(as.character(`...1`), "_"), `[`, 1)) %>% #  split '...1' and keep the first part
    pull(gene) # Extract the 'gene' column as a vector
  
  # Fetch chr
  chromosome_numbers <- sapply(genes, get_chromosome_number)
  gene_data <- data.frame(
    Gene = names(chromosome_numbers),
    Chromosome = chromosome_numbers,
    stringsAsFactors = FALSE)
  
  # Check for NA values in the 'CHR' column
  has_na_in_chr <- any(is.na(gene_data$Chromosome))
  # Print the result
  if (has_na_in_chr) {
    print("There are NA values in CHR column")
  } else {
    print("No NA values in CHR column")
  }
  
  # Count genes per chr
  genes_chr <- gene_data %>%
    group_by(Chromosome) %>%
    summarise(gene_count = n_distinct(Gene))
  
  # Plotting
  vg <- ggplot(genes_chr, aes(x = Chromosome, y = gene_count, fill = "#FF6666")) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    scale_x_discrete(limits = c(as.character(1:22), "X", "Y", "MT")) +
    scale_y_continuous(limits = c(0, 11), breaks = seq(0, 11, by = 1)) +
    labs(x = "Chromosome", y = "Number of Validated Genes") +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey", linewidth = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(size = 12),  
          axis.text.y = element_text(size = 12))
  
  # save  plot
  ggsave(plot = vg, filename = file.path(resultsDir, paste0("/functional/Internal_validation/barplot_chr_", name,".png")),  width = 10, height = 8)
}

# Create donut plots for each contrast
# Create Inner ring data frame
data_inner <- data.frame(
  labels = c("Validated", "Remaining Inner"),
  values = c(length(significant), (dim(unique_gene_positions)[1] - length(significant))), 
  color = c("#FF6666", "#FFB266"))
data_inner <- data.frame(
  labels = c("Validated", "Remaining Inner"),
  values = c(33, 210), 
  color = c("#FF6666", "#FFB266"))
# Outer ring data frame
dicovered_caas <- data %>% 
  distinct(Gene, Position) %>%
  nrow()

dicovered_genes <- data %>% 
  distinct(Gene, Position) %>%
  distinct(Gene) %>%
  nrow()

data_outer <- data.frame(
  labels = c("Detected", "Remaining Outer"),
  values = c(dicovered_genes, (16133 - dicovered_genes)),  # If total is 16133, remaining would be 16133 - dicovered_genes
  color = c("#66B266", "#FFFF66"))

# Define plotly figure
fig <- plot_ly() %>%
  add_pie(data = data_outer, labels = ~labels, values = ~values,
          domain = list(x = c(0, 1), y = c(0, 1)),
          name = "Outer", textinfo = 'value+percent', hole = 0.6,
          marker = list(colors = ~color), textfont = list(size = 20)) %>%  # Increase text size 
  add_pie(data = data_inner, labels = ~labels, values = ~values,
          domain = list(x = c(0.2, 0.8), y = c(0.2, 0.8)),  # Centered domain for the inner pie
          name = "Inner", textinfo = 'value+percent', hole = 0.5,  # Larger hole for the inner pie
          marker = list(colors = ~color), textfont = list(size = 15))

# Configure layout
fig <- fig %>% layout(title = '', showlegend = FALSE)

use_python("/home/vant/mambaforge/envs/tfm_r/bin/python3", required = TRUE)
save_image(fig,  file.path(resultsDir, "functional/Internal_validation/donuts_all4fam.png"), dpi = 300)
save_image(fig,  file.path(resultsDir, "functional/Internal_validation/donuts_cebi_atel.png"), dpi = 300)
save_image(fig,  file.path(resultsDir, "functional/Internal_validation/donuts_cebi_lemu.png"), dpi = 300)
save_image(fig,  file.path(resultsDir, "functional/Internal_validation/donuts_lemu_atel.png"), dpi = 300)
save_image(fig,  file.path(resultsDir, "functional/Internal_validation/donuts_cerco_cebi.png"), dpi = 300)
save_image(fig,  file.path(resultsDir, "functional/Internal_validation/donuts_cerco_atel.png"), dpi = 300)
save_image(fig,  file.path(resultsDir, "functional/Internal_validation/donuts_cerco_lemu.png"), dpi = 300)

# Import png image

image1 <- image_read(file.path(resultsDir,"functional/Internal_validation/donuts_cebi_lemu.png"))
image2 <- image_read(file.path(resultsDir,"functional/Internal_validation/cebi_lemu/boxplot_ARMT1_189.png"))
image3 <- image_read(file.path(resultsDir,"functional/Internal_validation/barplot_chr_cebi_lemu.png"))

# Convert png images to ggplot objects:
figA <- ggdraw() + draw_image(image1)
figB <- ggdraw() + draw_image(image2)
figC <- ggdraw() + draw_image(image3)

# Create the combined plot with figA on the left and the other two stacked on the right
# lemu cebi
combined_plot <- (figA | figB) / figC

# Define the layout to make figA larger and add the labels "A" and "B"
final_plot <- combined_plot  +
  plot_annotation(tag_levels = 'A')
final_plot

ggsave(final_plot,filename = file.path(resultsDir,"/functional/Internal_validation/internalval_cebi_lemu.png"),  dpi =  600)

# prepare supplementary figures
contrast <- list(all4fam = "CD164L2_140", cebi_atel = "CABLES1_174", cerco_atel = "HJURP_97", cerco_cebi = "DNAH8_2868", cerco_lemu = "MPP2_250", lemu_atel = "FREM2_2398")

# Loop
for (name in names(contrast)) {
  cat(paste(name), "contrast")
  image1 <- image_read(file.path(resultsDir, paste0("functional/Internal_validation/donuts_", name,".png")))
  image2 <- image_read(file.path(resultsDir, paste0("functional/Internal_validation/", name,"/boxplot_", contrast[[name]],".png")))
  
  # Convert png images to ggplot objects:
  figA <- ggdraw() + draw_image(image1)
  figB <- ggdraw() + draw_image(image2)
  
  final_plot <- (figA | figB)  +
    plot_layout(heights = c(1.5, 1)) +
    plot_annotation(tag_levels = 'A')
  
  ggsave(final_plot,filename = file.path(resultsDir,paste0("/functional/Internal_validation/internalval_", name,".png")), width = 6, height = 3, dpi =  400)
}
  
#---------------

# Load trasvar coordinates file
transvar_coord <- read_tsv(file.path(resultsDir,"/validation_outputs/coordinates_transvar.tab"), col_names = F)
# Load cebi_lemu list of internally validated genes
internal_valCAAS <- read_tsv(file.path(resultsDir,"functional/Internal_validation/cebi_lemu/cebi_lemu_genes.tab"), col_names = T)

cebi_lemu <- read_delim(file.path(resultsDir,"/functional/cebi_lemu_byFGN.txt"),col_names = TRUE ,
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)
# Give name to columns
colnames(transvar_coord) <- c("Gene",	"Trait",	"Position",	"transvar_coord")

merged_df <- inner_join(transvar_coord, internal_valCAAS, by = c("Gene", "Position"))  %>%
  dplyr::select(Gene,Position,transvar_coord)

merged_df <- left_join(merged_df, cebi_lemu, by = c("Gene", "Position"))  %>%
  dplyr::select(Gene,Position, Substitution,transvar_coord)

write.table(merged_df, file.path(resultsDir,"/functional/Internal_validation/cebi_lemu/internal_val_transvar_coord.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = T)
