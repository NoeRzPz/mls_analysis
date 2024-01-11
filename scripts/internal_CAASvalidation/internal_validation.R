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

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out/functional")

# Import df with CAAS and LQ values
all4fam <- read_tsv(file.path(resultsDir,"/4fam_sppclassification_results.tab"), col_names = T)
cebi_atel <- read_tsv(file.path(resultsDir,"/cebi_atel_sppclassification_results.tab"), col_names = T)
cebi_lemu <- read_tsv(file.path(resultsDir,"/cebi_lemu_sppclassification_results.tab"), col_names = T)
lemu_atel <- read_tsv(file.path(resultsDir,"/lemu_atel_sppclassification_results.tab"), col_names = T)
cerco_cebi <- read_tsv(file.path(resultsDir,"/cerco_cebi_sppclassification_results.tab"), col_names = T)
cerco_atel <- read_tsv(file.path(resultsDir,"/cerco_atel_sppclassification_results.tab"), col_names = T)
cerco_lemu <- read_tsv(file.path(resultsDir,"/cerco_lemu_sppclassification_results.tab"), col_names = T)

# Define path to python for saving plotly plots
use_python("~/mambaforge/envs/tfm_r/bin/python3", required = TRUE)

# Define functions 
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
# Define function to find outliers
findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

# Save 7 family contrast in a list called data
datas <- list(all4fam = all4fam, cebi_atel = cebi_atel, cebi_lemu = cebi_lemu, cerco_atel = cerco_atel, cerco_cebi = cerco_cebi, cerco_lemu = cerco_lemu, lemu_atel = lemu_atel)
 
# Loop
for (name in names(datas)) {
  cat(paste(name), "contrast\n")
  # Load phylogenetic tree 
  primates_tree <- read.tree(file.path(dataDir,"science.abn7829_data_s4.nex.tree"))
  
  # Find common names between tree labels and trait dataframe
  common_names <- intersect(primates_tree$tip.label, unique(datas[[name]]$label))
  
  data <- datas[[name]] %>%
    filter(label %in% common_names)
  
  # Filter the original data frame
  filtered_data <- data %>%
    group_by(Gene, Position) %>%
    # Filter genes with more than 5 cases of each level of the factor type_LQ
    filter(sum(type_LQ == "Long") >= 5, sum(type_LQ == "Short") >= 5)
  
  # change the reference level of a factor
  #gene_data$type_LQ <- relevel(factor(gene_data$type_LQ), ref = "Short")
  
  # Create unique gene-position pairs
  unique_gene_positions <- unique(filtered_data[c("Gene", "Position")])
  print(dim(unique_gene_positions)[1])
  
  
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
    
    gene_position_data <- gene_position_data %>%
      group_by(type_LQ) %>%
      mutate(outlier = ifelse(findoutlier(LQ), label, NA))
    
    # Make a name vector with the LQ values from or species
    #LQ <- gene_position_data$LQ
    #names(LQ) <- gene_position_data$label
    # Set a seed for reproducibility
    set.seed(123)
    
    result <- phylwilcoxtest(tree = pruned_tree, y = "LQ", x = "type_LQ", nsim = 3000, data = gene_position_data)
    
    # Store the results
    results[[paste(gene, position, sep = "_")]] <- result
    
    print(paste(paste(gene, position, sep = "_"), "processed"))
    
    pvalue <- signif(result$Pval, digits = 3)
    pvalue <- if (result$Pval < 0.001) {
      "< 0.001"
    } else {
      paste("=", round(result$Pval, digits = 3))
    }
    y.position <- (max(gene_position_data$LQ, na.rm = TRUE) * 1.05)
    
    # Plot of trait values by predictor categories
    LQplot <- ggplot(gene_position_data, aes(x = type_LQ, y = LQ, fill = type_LQ)) +
      geom_boxplot() +
      geom_text(aes(label = outlier), na.rm = TRUE, hjust = -.1, size = 5) +
      geom_point(size = 2, alpha = 0.5) +  # Add jitter to avoid overplotting
      scale_fill_manual(values = c("Short" = "#CC313D", "Long" = "#189ab4")) +
      labs(x = "CAAS Predictor", y = "LQ Value", title = paste( gene, " gene at position", position)) +
      # Draw the line between the two groups
      geom_segment(x = 1, xend = 2, y = y.position, yend = y.position, color = "black") + 
      annotate("text", x = 1.5, y = (y.position * 1.03), label = paste("p", pvalue), size = 7) +
      theme_classic() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.x = element_text(size = 18),  
            axis.title.y = element_text(size = 18),  
            plot.title = element_text(size = 20))
    
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
    write.table(cbind(names(pval), pval), file.path(resultsDir, paste0("/Internal_validation/",name,"/phylwilcoxon_pvalues.txt")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  }
  significant <- pval[which(pval <= 0.05)]
  #significant 
  # Adjust p-values using Benjamini-Hochberg method (not correct because the different positions in the same gene are not independent, think alternative)
  padj <- p.adjust(pval, method = "BH")
  #padj[which(padj <= 0.05)]
  # Adjust p-values for FDR (these are the q-values)
  q_values <- p.adjust(pval, method = "fdr")
  # Display q-values
  #q_values[which(q_values <= 0.05)]
  
  # loop through the list to print significant plots
  for (gene_pos in names(significant)) {
    print(graphs[[gene_pos]])
    # save all significant plots to files
    ggsave(filename = file.path(resultsDir,paste0("/Internal_validation/",name,"/boxplot_", gene_pos, ".png")), plot = graphs[[gene_pos]], width = 10, height = 8)
  }
  # Create donut plots for each contrast
  # Create Inner ring data frame
  data_inner <- data.frame(
    labels = c("Validated", "Remaining Inner"),
    values = c(length(significant), (dim(unique_gene_positions)[1] - length(significant))), 
    color = c("#FF6666", "#FFB266"))
  
  # Outer ring data frame
  dicovered_caas <- datas[[name]] %>% 
    distinct(Gene, Position) %>%
    nrow()
  
  dicovered_genes <- datas[[name]] %>% 
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
            marker = list(colors = ~color), textfont = list(size = 17), textposition = 'inside')
  
  # Configure layout
  fig <- fig %>% layout(title = '', showlegend = FALSE)
  
  
  save_image(fig,  file.path(resultsDir, paste0("/Internal_validation/donuts_",name,".png")), dpi = 300, width = 2, height = 3)
}

#hist(pval, breaks = 50)

#---------------
# Create barplots to spot chromosomes where internally validated genes fall

# List of internaly validated genes you want to query
all4fam_validated <- read_delim(file.path(resultsDir,"/Internal_validation/all4fam/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cebi_atel_validated <- read_delim(file.path(resultsDir,"/Internal_validation/cebi_atel/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cebi_lemu_validated <- read_delim(file.path(resultsDir,"/Internal_validation/cebi_lemu/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cerco_atel_validated <- read_delim(file.path(resultsDir,"/Internal_validation/cerco_atel/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cerco_cebi_validated <- read_delim(file.path(resultsDir,"/Internal_validation/cerco_cebi/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
cerco_lemu_validated <- read_delim(file.path(resultsDir,"/Internal_validation/cerco_lemu/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)
lemu_atel_validated <- read_delim(file.path(resultsDir,"/Internal_validation/lemu_atel/phylwilcoxon_pvalues.txt"), delim = "\t",   col_names = TRUE)

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
  cat(paste(name), "contrast\n")
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
          axis.text.x = element_text(size = 16),  
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))
  
  # save  plot
  ggsave(plot = vg, filename = file.path(resultsDir, paste0("/Internal_validation/barplot_chr_", name,".png")), dpi =  300, width = 10, height = 8)
}


# Import png image
image1 <- image_read(file.path(resultsDir,"/Internal_validation/donuts_cebi_lemu.png"))
image2 <- image_read(file.path(resultsDir,"/Internal_validation/cebi_lemu/boxplot_CUBN_2476.png"))
image3 <- image_read(file.path(resultsDir,"/Internal_validation/barplot_chr_cebi_lemu.png"))

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

ggsave(final_plot,filename = file.path(resultsDir,"/Internal_validation/internalval_cebi_lemu.png"), dpi =  600, width = 7, height = 8)

# prepare supplementary figures
contrast <- list(all4fam = "CD164L2_140", cebi_atel = "CABLES1_174", cerco_atel = "HJURP_97", cerco_cebi = "DNAH8_2868", cerco_lemu = "MPP2_250", lemu_atel = "FREM2_2398")

# Loop
for (name in names(contrast)) {
  cat(paste(name), "contrast\n")
  image1 <- image_read(file.path(resultsDir, paste0("/Internal_validation/donuts_", name,".png")))
  image2 <- image_read(file.path(resultsDir, paste0("/Internal_validation/", name,"/boxplot_", contrast[[name]],".png")))
  
  # Convert png images to ggplot objects:
  figA <- ggdraw() + draw_image(image1)
  figB <- ggdraw() + draw_image(image2)
  
  final_plot <- (figA | figB)  +
    plot_layout(heights = c(1.5, 1)) +
    plot_annotation(tag_levels = 'A')
  
  ggsave(final_plot,filename = file.path(resultsDir,paste0("/Internal_validation/internalval_", name,".png")), width = 5, height = 3, dpi =  400)
}
  
