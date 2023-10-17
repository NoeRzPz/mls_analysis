# Clean environment
rm(list = ls())

# Load packages
library(tidyverse, quietly = TRUE)
library(taxize)

# Directories setting up
workingDir <- getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir,"out") 

# Load info
spp_names <- read_csv(file.path(resultsDir,"spp_names.csv"), col_names = TRUE)

# Retrieve taxonomic family using tax_name
spp_family <- tax_name(spp_names$x, get = "family",  db = "ncbi") #accepted = TRUE,

# Check if there are missing families
na_fam <- sum(is.na(spp_family$family))
if (na_fam > 0) {
  spp_family <- spp_family %>%
    mutate( genus = word(query, 1, 1),
            family = ifelse(is.na(family), first(family[!is.na(family)]), family)) %>%
    ungroup()
  cat("NA family has been replaced by the equivalent genus with family\n")
} else {
  cat("No missing families in dataset\n")
}

spp_family <- spp_family %>%
  rename(Family = family) 

# Save info obtained
write.csv(spp_family, file.path(dataDir,"spp_family.csv"), row.names = FALSE)
