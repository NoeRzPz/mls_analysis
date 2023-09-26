# Functional enrichment detection methods to find functionally-related genomic elements experiencing convergent evolutionary rates as a group

# Load package
library(RERconverge)
library(ggplot2)
library(ggrepel)
# Import inputs: 
# output from RERconverge correlation functions 
res1 <- read.csv("out/continuostrait_cor_res.csv", row.names = 1) 
# Import Pathway Annotations
annots <- read.gmt("data/c2.all.v6.2.symbols.gmt")

# RERconverge enrichment functions expect pathways to be in named pathway-group lists contained within a list 
annotlist <- list(annots)
names(annotlist) <- "MSigDBpathways"

# Create a ranked gene list
stats <- getStat(res1)

# Calculate Enrichment
# Use Wilcoxon Rank-Sum Test on a list of genes ranked based on their correlation statistics
enrichment <- fastwilcoxGMTall(stats, annotlist, outputGeneVals = T, num.g = 10)

# Visualize results
# Sort the data by AdjustedPValue for better visualization
signif_res <- enrichment$MSigDBpathways %>%
  arrange(pval) %>%
  filter(pval <= 0.05) 

# Create a bar plot
ggplot(signif_res, aes(x = reorder(row.names(signif_res), p.adj), y = stat, fill = p.adj)) +
  geom_bar(stat = "identity") +
  labs(x = "Gene Set", y = "Enrichment Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient(low = "blue", high = "red") +
  ggtitle("Enrichment Analysis Results") +
  scale_y_continuous(expand = c(0, 0))

# Create a scatterplot
top_terms <- head(signif_res, 5)
ggplot(enrichment$MSigDBpathways, aes(x = stat, y =  -log10(p.adj))) +
  geom_point(aes(size = num.genes), alpha = 0.5) +
  geom_text_repel(data = top_terms, aes(label = row.names(top_terms)), 
                  box.padding = 0.5, point.padding = 0.2, 
                  max.overlaps = Inf, nudge_x = 0.1, nudge_y = 0.1,
                  size = 3) +
  labs(x = "Statistic (stat)", y = "-log10(Adjusted P-Value)") +
  scale_size_continuous(range = c(2, 10)) +  # Adjust the size range as needed
  theme_minimal() +
  ggtitle("Enrichment Analysis Results")

# Tras este script quedaria: use branch-site models to determine if fast-evolving genes are experiencing relaxation of constraint or undergoing directional selection