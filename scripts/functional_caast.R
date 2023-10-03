# Clean environment
rm(list = ls())

# Load packages
library(readr)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler, quietly = TRUE) 

# Load data 
# I start with the discovery dataframe (I have to use data after bootstraping to filter them )
longevity_farre.discov <- read_delim("data/longevity.farre2021.nofilter.tab", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
longevity_farre.discov$adj.pval <- p.adjust(longevity_farre.discov$Pvalue, method = "BH",  n = length(longevity_farre.discov$Pvalue))
longevity_sigf.05 <- longevity_farre.discov %>%
  filter(adj.pval <= 0.05)

head(longevity_sigf.05)

# Query to Translate RefSeq to Ensbl, Entrez:
refseq_id <- unique(longevity_sigf.05$Gene)

# With bitr I obtain desired anotations
ann <- bitr(refseq_id, fromType = "REFSEQ", toType =  c("REFSEQ","ENSEMBL","ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
background_genes <- bitr(unique(longevity_farre.discov$Gene), fromType = "REFSEQ", toType =  c("REFSEQ","ENSEMBL","ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
#the biomaRt package queries the Ensembl database, which might not always have the latest NCBI data. always cross-check your results and update the BioMart data when needed.

# OVER REPRESENTATION ANALYSIS (ORA)
  
#  To determine what a priori defined gene sets  are more represented that we could expect by chance in our subset of significant genes 

# We provide ENSEMBL IDs of our significant genes.  Specify option readable TRUE to translate ENSEMBL  IDs to gene symbols

#Especificamos que calcule el ORA para los  procesos biológicos (BP) de GO, e empleamos minGSSize y maxGSSize para restringir el tamaño de los gene sets.
#minGSSize, es el tamaño mínimo de genes anotados a un *Ontology term* que será testado (nos interesa decartar los muy pequeños porque alcanzarían fácilmente significancia). maxGSSize es el tamaño máximo de genes anotados que será testado (nos interesa descartar los muy grandes porque serán demasiado generales).

# Perform left join
#longevity_sigf.05 <- left_join(longevity_sigf.05, ann , by = c("Gene" = "REFSEQ"), multiple = "all")

set.seed(123)  #fijamos la semilla de aleatorización, evitamos que cada vez el resultado pueda variar

ego1 <- enrichGO(gene          = ann$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 minGSSize = 15,
                 maxGSSize = 500,
                 pvalueCutoff  = 0.05, #0.01
                 qvalueCutoff  = 0.05,
                 readable = TRUE, 
                 universe = background_genes$ENSEMBL)

dim(ego1) # con dim podemos saber el número de BP significativos 


#Empleamos el método `simplify` para eliminar términos redundantes. con 'select_fun' seleccionamos un término representativo de los términos reduntantes por su menor pvalor ajustado. Los términos redundantes serán los que tengan una similitud mayor que 0.7. 

ego1 <- simplify(ego1 , cutoff = 0.7, by = "p.adjust", select_fun = min)

# Indica el número de BP GOs  significativos
dim(ego1)

# number of BP significant
sig.BP <- head(ego1, nrow(ego1)) 
head(sig.BP)  

#Visualizamos los resultados de forma gráfica, para tener una mejor visión global de los principales procesos afectados por estos genes
#`goplot` nos devuelve un DAG de los términos GO significativos.

goplot(ego1)

dotplot(ego1, showCategory = 20)

ego1p <- enrichplot::pairwise_termsim(ego1)
emapplot(ego1p, cex_label_category = 0.8, showCategory = 20)


cnetplot(ego1, showCategory = 15, cex_label_category = 0.8)


# KEGG pathway over-representation analysis
kk <- enrichKEGG(gene         = ann2$ENSEMBL,
                 organism     = 'hsa',
                 universe = as.character(background_genes$ENSEMBL),
                 pvalueCutoff = 0.05)
# Ideas: 
# Study if those genes share a specific regulation eg:if they have binding sites to common transcription factors
# Study if they are regulated by the same miRNAs. 

sessionInfo()
