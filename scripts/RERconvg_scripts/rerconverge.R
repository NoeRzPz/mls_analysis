# Clean environment
rm(list = ls())

# Loading packages
library(readr)
library(tidyverse)
library(RColorBrewer)
library(RERconverge)
library(ape)
library(phytools)

# Setting up directories
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Feeding data
#genes_tree <- readTrees(file.path(dataDir,"/ALL_FEB23_geneTrees.txt"))
#saveRDS(genes_tree, file = "out/gene_trees.rds") 
gene_trees <- readRDS(file = file.path(resultsDir,"/gene_trees.rds"))
primates_tree <- read.tree(file.path(dataDir,"/caas_tool_inputs/primates.tree.nwk"))
primates_LQ <- read_delim(file.path(resultsDir,"/descriptive/trait_LQ.tab"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Make a name vector with the LQ values from or species
primates_lq <- primates_LQ$LQ
names(primates_lq) <- primates_LQ$label

# Make a name vector with the LQ values from or species
#primates_mls <- primates_MLSresiduals$pgls_residuals
#names(primates_mls) <- primates_MLSresiduals$label

# We represent the continuous trait in the primates phylogeny

# Find the common names between vector names and the names in the trees 
common_names <- intersect(gene_trees$masterTree$tip.label, names(primates_lq))
# Drop the tips that are NOT in common_names
sub_tree <- drop.tip(gene_trees$masterTree, setdiff(gene_trees$masterTree$tip.label, common_names))

# Phylogeny trait display with phytools
dotTree(sub_tree, primates_lq,length = 10,ftype = "i")

# Another one to display the trait
setdiff(names(primates_lq), sub_tree$tip.label)
primates_LQ <- primates_lq[sub_tree$tip.label]


# Define a color palette function for a single color with varying intensities
# Adjust the number of colors (n) to control intensity levels
palette <- brewer.pal(n = 9, name = "Blues")

# Create a continuous palette with colorRampPalette
color_palette <- colorRampPalette(palette)
colors <- color_palette(100)

# Create the contMap object
obj <- contMap(sub_tree, primates_LQ, plot = FALSE)
obj <- setMap(obj, colors = colors)

png(filename = file.path(resultsDir,"/descriptive/LQ_phylogeny.png"),  width = 20, height = 16, units = "cm", res = 300)

#plot(obj, lwd = 3, xlim = c(-0.2,3.6))
plot(obj,legend = FALSE, lwd = 4, ylim = c(1-0.09*(Ntip(obj$tree)-1), Ntip(obj$tree)), mar = c(5.1,0.4,0.4,0.4))
add.color.bar(0.5,obj$cols,title = "LQ Tree Distribution",
              lims = obj$lims, digits = 3, prompt = FALSE, x = 0,
              y = 1-0.08*(Ntip(obj$tree)-1), lwd = 4, fsize = 1,subtitle = "")
# add  x-axis
axis(1)
title(xlab = "Time from the root")

dev.off()

# Calculating RERs and plot results to check for heteroscedasticity correction
# Create and save plots as PNG
#png("out/heterocedasticity_correction.png", width = 800, height = 600, res = 300)  
mamRERw <- getAllResiduals(gene_trees, useSpecies = names(primates_lq), transform = "sqrt", weighted = T, scale = T, plot = TRUE)
#, cutoff = 0 cutoff 0 for no cutoff
#mamRERw <- getAllResiduals(gene_trees, transform = "sqrt", weighted = T, scale = T, plot = TRUE)
#dev.off()  # Close the PNG device to save the file

#Save object for later use
#saveRDS(mamRERw, file = file.path(resultsDir,"/mamRERw.rds"))
mamRERw <- readRDS(file = file.path(resultsDir,"/mamRERw.rds"))

# Continuous Trait Analysis
###########################

#mamRERwMLS <- getAllResiduals(gene_trees, useSpecies = names(primates_mls), transform = "sqrt", weighted = T, scale = T, plot = TRUE, cutoff = 0)
#saveRDS(mamRERwMLS, file = file.path(resultsDir,"/mamRERwMLS.rds"))
#charpaths <- char2Paths(primates_mls, gene_trees)
#res <- correlateWithContinuousPhenotype(mamRERwMLS, charpaths, min.sp = 10, winsorizeRER = 1, winsorizetrait = 1)#ni con 3 ni con 2 sale signnificativos


# Convert the trait vector to paths comparable to the paths in the RER matrix
charpaths <- char2Paths(primates_lq, gene_trees)
# metric diff is used, so branch lengths are assigned to the trait tree based on the difference in trait values on the nodes connected to that branch
#tree2Paths
#to find correlations between the rate of evolution of a genomic element (encoded in the RER matrix) and the rate of change of a phenotype (encoded in charpaths). The final output is the list of input genes with relevant statistics. As input, I provide the RER matrix and trait path.
#I set method=“p” to use a Pearson correlation, which is the most appropriate option for continuous traits. I set min.pos=0 (min.pos for “minimum positive”) so no foreground species are required; correlation statistics will be calculated even if all species have negative RER values. Finally, winsorize=3 pulls the outermost three points in all four directions to the center of the data before calculating correlations to mitigate the effects of outlier points.
res <- correlateWithContinuousPhenotype(mamRERw, charpaths, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
write.csv(res, file.path(resultsDir,"rerconverge/continuostrait_cor_res_standardfilters.csv"), row.names = TRUE)

head(res[order(res$P),])

# Antes lo corri con min.sp = 10, winsorizeRER = 3, winsorizetrait = 3, no salia nada y por eso probé relajando las opciones con uno si salen, por homo sapiens
#res <- correlateWithContinuousPhenotype(mamRERw, charpaths,min.sp = 10, winsorizeRER = 1, winsorizetrait = 1)

# Histogram pvalues
png(file.path(resultsDir,"rerconverge/histogram_pval.png"), width = 1500, height = 1200, res = 300)
hist(res$P, breaks = 15,  main = "P-values for correlation between genes and LQ", xlab = "Pvalues")
dev.off()

# Order by pval & Filter genes with adjusted p-value <= 0.05
signif_res_0.05 <- res %>%
  arrange(P) %>%
  filter(P <= 0.05) 
write.csv(signif_res_0.05, file.path(resultsDir,"rerconverge/signif_lq_genes_005.csv"), row.names = TRUE)


# Plotting RER for gene 'SLC25A6' using the sorted data
# for (i in rownames(signif_res_0.05)) {
#   png(file.path(resultsDir,paste0("/rerconverge/", i, "_RER_plot.png")), width = 9, height = 7, units = "in", res = 300)
#   plotRers(mamRERw, i , phenv = primates_lq)
#   dev.off()
# }


# Check the impact of one or a few species on the overall correlation. To assess this risk, we examine individual gene correlation plots:
for (i in rownames(signif_res_0.05)) {
  png(file.path(resultsDir,paste0("/rerconverge/", i, "_correlation_plot.png")), width = 9, height = 7, units = "in", res = 300)
  x <- charpaths
  y <- mamRERw[i,]
  pathnames <- namePathsWSpecies(gene_trees$masterTree) 
  names(y) <- pathnames
  plot(x,y, cex.axis = 1, cex.lab = 1, cex.main = 1, xlab = "LQ Change", 
       ylab = "Evolutionary Rate", main = paste("Gene", i, "Pearson Correlation"),
       pch = 19, cex = 1, xlim = c(-2,2))
  text(x,y, labels = names(y), pos = 4)
  abline(lm(y~x), col = 'red',lwd = 3)
  dev.off()
}
