# Loading packages
if (!require("RERconverge", character.only = T, quietly = T)) {
  require(devtools)
  install_github("nclark-lab/RERconverge", ref = "master") 
  #"ref" can be modified to specify a particular branch
}
library(RERconverge)
#library(ggplot2)

# Setting up files directory
rerpath = find.package('RERconverge') 
print(rerpath)

# Feeding data
treefile = "subsetMammalGeneTrees.txt" 
Trees <- readTrees(paste(rerpath,"/extdata/",treefile,sep = ""))

data("logAdultWeightcm")

# Calculating RERs and plot results to check for heteroscedasticity correction

# Create and save plots as PNG
png("out/heterocedasticity_correction.png", width = 800, height = 600, res = 300)  
mamRERw <- getAllResiduals(Trees,useSpecies = names(logAdultWeightcm), transform = "sqrt", weighted = T, scale = T, plot = TRUE)
dev.off()  # Close the PNG device to save the file

#Save object for later use
saveRDS(mamRERw, file = "out/mamRERw.rds") 
#newmamRERw = readRDS("out/mamRERw.rds")


# Binary Trait Analysis
#######################
noneutherians <- c("Platypus","Wallaby","Tasmanian_devil","Opossum")

# Predefined foregrounds in a file
marineb <- read.tree(paste(rerpath,"/extdata/MarineTreeBinCommonNames_noCGM.txt",sep = ""))
marinebrooted = root(marineb,outgroup = noneutherians)
mb1 = marineb
mb1$edge.length = c(rep(1,length(mb1$edge.length)))

# Representing the tree
binplot1 <- plotTreeHighlightBranches(mb1, outgroup = noneutherians,
                                   hlspecies = which(marineb$edge.length == 1), hlcols = "blue", 
                                   main = "Foreground branches highlighted") 
# Make average and gene tree plots
# Plot average tree
# Save tree plots
png("out/average_tree.png", width = 800, height = 600)
# Branch lengths on tree represent average rates across all genes
avgtree <- plotTreeHighlightBranches(Trees$masterTree, outgroup = noneutherians,
                                     hlspecies = marineextantforeground, hlcols = "blue", 
                                     main = "Average tree") 
dev.off()

# Other options to define foreground
marineb2a = foreground2Tree(marineextantforeground, Trees, clade = "ancestral", useSpecies = names(logAdultWeightcm))
marineb2b = foreground2Tree(marineextantforeground, Trees, clade = "terminal", useSpecies = names(logAdultWeightcm))
marineb2c = foreground2Tree(marineextantforeground, Trees, clade = "all", useSpecies = names(logAdultWeightcm))
marineb2d = foreground2Tree(marineextantforeground, Trees, clade = "all", weighted = TRUE, useSpecies = names(logAdultWeightcm))

#plot(marineb, main = "Manually specified binary tree")
# If we want to select foreground branches by hand:
#marineb3 <- click_select_foreground_branches(Trees$masterTree)
#marineb3rooted <- root(marineb3,outgroup=c("Platypus", "Wallaby","Tasmanian_devil","Opossum"))

# For genes may not have data for all species, meaning that their phylogeny will be a subset of the full phylogeny. To plot RERs for these genes and to correlate them with trait evolution, 
phenvMarine <- tree2Paths(marineb, Trees)

marineextantforeground = c("Walrus","Seal","Killer_whale","Dolphin","Manatee")
phenvMarine2 <- foreground2Paths(marineextantforeground, Trees, clade = "all") # "terminal" can be used instead too

phenvMarine2b <- tree2Paths(marineb2b, Trees)

# Correlate gene evolution with binary trait evolution
corMarine <- correlateWithBinaryPhenotype(mamRERw, phenvMarine, min.sp = 10, min.pos = 2, weighted = "auto")

# Overall pattern of association across all genes
hist(corMarine$P, breaks = 15, xlab = "Kendall P-value",  main = "P-values for correlation between genes and marine environment")

# Order by pval & Filter genes with adjusted p-value <= 0.05
signif_res_0.05 <- corMarine %>%
  arrange(P) %>%
  filter(p.adj <= 0.05) 
write.csv(signif_res_0.05, "out/signif_longevity_genes_005.csv", row.names = TRUE)

# Order by pval & Filter genes with adjusted p-value <= 0.1
signif_res_0.1 <- corMarine %>%
  arrange(P) %>%
  filter(p.adj <= 0.1) 
write.csv(signif_res_0.1, "out/signif_longevity_genes_010.csv", row.names = TRUE)

# Visualize &  save RERs for significant genes (plot RERs for each gene):
for (i in rownames(signif_res_0.1)) {
  png(paste0("out/RER_", i, ".png"), width = 7, height = 5, units = "in", res = 300)
  plotRers(mamRERw, i ,phenv = phenvMarine)
  dev.off()
}
#visualize RERs along branches as a heatmap
for (i in rownames(signif_res_0.1)) {
  png(paste0("out/heat_RER_", i, ".png"), width = 9, height = 7, units = "in", res = 300)
  newrers = treePlotRers(treesObj = Trees, rermat = mamRERw, index = i, type = "c", nlevels = 9, figwid = 10)
  dev.off()
}

# Save individual significant gene trees
# Branch lengths on tree represent rates specifically for the gene
for (i in rownames(signif_res_0.1)) {
  png(paste0("out/", i, "_tree.png"), width = 7, height = 8, units = "in", res = 300)
  genetree <- plotTreeHighlightBranches(Trees$trees[[i]], outgroup = noneutherians, 
                            hlspecies = marineextantforeground, hlcols = "blue", 
                            main = paste(i, "tree")) # which(phenvMarine == 1)
  dev.off()
}

# Save & plot RERs as tree for that gene
for (i in rownames(signif_res_0.1)) {
  png(paste0("out/tree_RER_", i, ".png"), width = 9, height = 7, units = "in", res = 300)
  gene_rers = returnRersAsTree(Trees, mamRERw, i, plot = TRUE, phenv = phenvMarine) 
  dev.off()
}
# Save RERs as tree for all genes
multirers = returnRersAsTreesAll(Trees, mamRERw)
write.tree(multirers, file = 'out/RERs_trees.nwk', tree.names = TRUE)

# Continuous Trait Analysis
###########################
# Vector names must match the names in the trees read in previously
data("logAdultWeightcm")
head(logAdultWeightcm)

# Convert the trait vector to paths comparable to the paths in the RER matrix
charpaths <- char2Paths(logAdultWeightcm, Trees)
# metric diff is used, so branch lengths are assigned to the trait tree based on the difference in trait values on the nodes connected to that branch
res <- correlateWithContinuousPhenotype(mamRERw, charpaths, min.sp = 10, winsorizeRER = 3, winsorizetrait = 3)
write.csv(res, "out/continuostrait_cor_res.csv", row.names = TRUE)

# Overall pattern of association across all genes
hist(res$P, breaks = 15)

# Order by pval & Filter genes with adjusted p-value <= 0.05
signif_res_0.05 <- res %>%
  arrange(P) %>%
  filter(p.adj <= 0.05) 
head(signif_res_0.05)
write.csv(signif_res_0.05, "out/signif_weight_genes_005.csv", row.names = TRUE)

# Order by pval & Filter genes with adjusted p-value <= 0.1
signif_res_0.1 <- res %>%
  arrange(P) %>%
  filter(p.adj <= 0.1) 
write.csv(signif_res_0.1, "out/signif_weight_genes_010.csv", row.names = TRUE)

# Check the impact of one or a few species on the overall correlation. To assess this risk, we examine individual gene correlation plots:
for (i in rownames(signif_res_0.1)) {
  png(paste0("out/", i, "_correlation_plot.png"), width = 9, height = 7, units = "in", res = 300)
  x <- charpaths
  y <- mamRERw[i,]
  pathnames <- namePathsWSpecies(Trees$masterTree) 
  names(y) <- pathnames
  plot(x,y, cex.axis = 1, cex.lab = 1, cex.main = 1, xlab = "Weight Change", 
     ylab = "Evolutionary Rate", main = paste("Gene", i, "Pearson Correlation"),
     pch = 19, cex = 1, xlim = c(-2,2))
  text(x,y, labels = names(y), pos = 4)
  abline(lm(y~x), col = 'red',lwd = 3)
  dev.off()
}


