library(data.table)
library(tidyverse)
library(ggplot2)

# Load v4 and v9 tables
# ASVs with taxonomy with confidence > 50% by idtaxa; only sequences seen in at least two different samples with at least 3 copies in the whole dataset

# original table from https://zenodo.org/records/13881418
v9 <- fread(file="datasets/TARA-Oceans_18S-V9_dada2_table.tsv.gz", sep="\t", quote=F)
v9 <- v9[grep(v9$taxonomy, pattern="^Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta"),]
v9$method <- 'V9'

# original table from https://zenodo.org/records/13881376
v4 <- fread(file="datasets/TARA-Oceans_18S-V4_dada2_table.tsv.gz", sep="\t", quote=F)
v4 <- v4[grep(v4$taxonomy, pattern="^Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta"),]
v4$method <- 'V4'

# Define a function to process the data
process_data <- function(data, samples) {
  # Select only relevant samples and compute read abundance
  samples <- samples[samples$sample_id_pangaea %in% colnames(data),]
  combined_data <- data %>%
    select(amplicon, taxonomy, confidence, total, spread, sequence, samples$sample_id_pangaea)
  
  combined_data$total_subset <- rowSums(combined_data %>% select(samples$sample_id_pangaea))
  
  # diatom ASVs with at least 3 counts
  #  combined_data <- combined_data[combined_data$total_subset != 0,]
  combined_data <- combined_data[combined_data$total_subset >= 3,]
  combined_data$total_subset <- NULL
  
  # diatom ASVs with ocurrance in at least 2 different samples
  sample_columns <- combined_data[, 7:ncol(combined_data)] #Identify the sample columns (7th to the last)
  occurrences <- rowSums(sample_columns > 0) #Calculate the number of samples where each species occurs
  rm(sample_columns)
  combined_data <- combined_data[occurrences >= 2, ] #Filter rows with occurrences in at least two different samples
  rm(occurrences)
  
  return(combined_data)
}

# Read sample metadata
# original table from https://zenodo.org/records/7229815
mysamples <- read.csv("datasets/context_general.tsv", sep = "\t")
mysamples <- mysamples %>%
  filter(depth %in% c("SRF", "DCM")) %>%
  filter(!size_fraction %in% c("0.22-1.6", "0.22-3", "0.8-3"))

# Process v4 and v9 data
v4_processed <- process_data(v4, mysamples)
v9_processed <- process_data(v9, mysamples)

#############################
# rarefaction curve #########
#############################

library('vegan')
v4_matrix <- v4_processed[,!c(1:6)]
v4_matrix <- t(v4_matrix) # species in columns, samples in rows
# Calculate the rarefaction curve using specaccum
rarefaction_model_V4 <- specaccum(v4_matrix, method="random") #method = "rarefaction"

v9_matrix <- v9_processed[,!c(1:6)]
v9_matrix <- t(v9_matrix) # species in columns, samples in rows
# Calculate the rarefaction curve using specaccum
rarefaction_model_V9 <- specaccum(v9_matrix, method="random")  #method = "rarefaction"

# Plot the first rarefaction curve
plot(rarefaction_model_V4, xlab="Number of Samples", ylab="Diatom ASV number", main="V4 & V9", xlim=c(1, max(nrow(v4_matrix), nrow(v9_matrix))), ylim=c(0, max(rarefaction_model_V4$richness, rarefaction_model_V9$richness)), col="blue", ci.col="gray")

# Add the second rarefaction curve
lines(rarefaction_model_V9, col="blue", ci.col="gray")

########
# estimate total abundance
specpool(v9_matrix)
specpool(v4_matrix)
#############################
#############################
