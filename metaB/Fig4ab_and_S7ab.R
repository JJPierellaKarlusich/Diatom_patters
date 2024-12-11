
library(data.table)
library(tidyverse)
library(ggplot2)
library(vegan)
library(cowplot)

#### Formating the data

# Load v4 and v9 tables
# ASVs with taxonomy with confidence > 50% by idtaxa; only sequences seen in at least two different samples with at least 3 copies in the whole dataset

# original table from https://zenodo.org/records/13881418
v9 <- fread(file="datasets/TARA-Oceans_18S-V9_dada2_table.tsv.gz", sep="\t", quote=F)
v9 <- v9[grep(v9$taxonomy, pattern="^Root;Eukaryota;"),]
v9$method <- 'V9'

# original table from https://zenodo.org/records/13881376
v4 <- fread(file="datasets/TARA-Oceans_18S-V4_dada2_table.tsv.gz", sep="\t", quote=F)
v4 <- v4[grep(v4$taxonomy, pattern="^Root;Eukaryota;"),]
v4$method <- 'V4'

# Load sample metadata
# original table from https://zenodo.org/records/7229815
mysamples <- read.csv("datasets/context_general.tsv", sep = "\t")
mysamples <- mysamples %>%
  filter(depth %in% c("SRF", "DCM")) %>%
  filter(!size_fraction %in% c("0.22-1.6", "0.22-3", "0.8-3"))


# Define a function to process the data
process_diatom_data <- function(data, samples) {
  
  # Filter diatoms
  combined_data <- data[grep(x=data$taxonomy, pattern = ";Bacillariophyta"),]
  
  # classify by diatom class
  
  #"Mediophyceae" # polar and 'non-polar' centric
  #"Coscinodiscophyceae" #Radial-centric-basal
  combined_data$taxonomy <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Radial-centric-basal-Coscinodiscophyceae.*', replacement = 'Coscinodiscophyceae')
  combined_data$taxonomy <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Raphid-pennate.*', replacement = 'Raphid_pennate')
  combined_data$taxonomy <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Polar-centric-Mediophyceae.*', replacement = 'Mediophyceae')
  combined_data$taxonomy <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Araphid-pennate.*', replacement = 'Araphid_pennate')
  
  combined_data$taxonomy <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Bacillariophyta_XX;Bacillariophyta_XXX;Bacillariophyta_XXX_sp.', replacement = 'unknown_class')
  combined_data$taxonomy <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;unclassified_Bacillariophyta_X', replacement = 'unknown_class')
  
  # Select only relevant samples 
  samples <- samples[samples$sample_id_pangaea %in% colnames(combined_data),]
  combined_data <- combined_data %>%
    select(amplicon, taxonomy, confidence, total, spread, sequence, samples$sample_id_pangaea)
  
  combined_data$total_subset <- rowSums(combined_data %>% select(samples$sample_id_pangaea))
  # ASVs with at least 3 counts in the whole selected dataset
  combined_data$total_subset <- rowSums(combined_data %>% select(samples$sample_id_pangaea))
  #combined_data <- combined_data[combined_data$total_subset != 0,]
  combined_data <- combined_data[combined_data$total_subset >= 3,]
  combined_data$total_subset <- NULL
  # ASVs with ocurrance in at least 2 different samples the whole selected 
  sample_columns <- combined_data[, 7:ncol(combined_data)] #Identify the sample columns (7th to the last)
  occurrences <- rowSums(sample_columns > 0) #Calculate the number of samples where each species occurs
  rm(sample_columns)
  combined_data <- combined_data[occurrences >= 2, ] #Filter rows with occurrences in at least two different samples
  rm(occurrences)
  
  # abundance sum per diatom class
  myabundance <- combined_data %>% select(taxonomy, samples$sample_id_pangaea)
  myabundance <- aggregate(data=myabundance, FUN=sum, . ~ taxonomy)
  row.names(myabundance) <- myabundance$taxonomy
  myabundance$taxonomy <- NULL
  myabundance <- t(myabundance)
  
  # Calculate total eukaryotic read abundance and normalize diatom class relative abundance
  euk_data <- data %>% select(samples$sample_id_pangaea)
  euk_data <- t(euk_data)
  ReadAbundance <- rowSums(euk_data, na.rm = TRUE)
  rm(euk_data)
  ReadAbundance <- data.frame(ReadAbundance_Euks = ReadAbundance)
  Indices <- merge(x=myabundance, y=ReadAbundance, by="row.names", all=T)
  rm(myabundance, ReadAbundance)
  colnames(Indices)[1] <- 'sample'
  Indices$Araphid_pennate <- 100 * Indices$Araphid_pennate / Indices$ReadAbundance_Euks
  Indices$Raphid_pennate <- 100 * Indices$Raphid_pennate / Indices$ReadAbundance_Euks
  Indices$Coscinodiscophyceae <- 100 * Indices$Coscinodiscophyceae / Indices$ReadAbundance_Euks
  Indices$Mediophyceae <- 100 * Indices$Mediophyceae / Indices$ReadAbundance_Euks
  Indices$unknown_class <- 100 * Indices$unknown_class / Indices$ReadAbundance_Euks
  
  # Merge with sample metadata
  Indices <- merge(Indices, samples, by.x="sample", by.y="sample_id_pangaea", all.x=TRUE)
  
  # calculate mean of technical replicates
  Indices <- aggregate(data=Indices, FUN=mean, cbind(ReadAbundance_Euks, Araphid_pennate, Raphid_pennate, Coscinodiscophyceae,Mediophyceae, unknown_class) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction)
  
  # Organize size fractions and mean equivalent size fractions
  Indices$size_fraction <- recode(Indices$size_fraction, "3-20" = "3/5-20", "5-20" = "3/5-20", "0.8-5" = "0.8-5/2000", "0.8->" = "0.8-5/2000", "0.8-20" = "0.8-5/2000")
  Indices$size_fraction <- factor(Indices$size_fraction, levels = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
  Indices <- aggregate(data=Indices, FUN=mean,  cbind(ReadAbundance_Euks, Araphid_pennate, Raphid_pennate, Coscinodiscophyceae,Mediophyceae, unknown_class) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction)
  
}



# Process v4 and v9 data
v4_processed <- process_diatom_data(v4, mysamples)
v4_processed$method <- 'V4'
v9_processed <- process_diatom_data(v9, mysamples)
v9_processed$method <- 'V9'
myabundance <- rbind(v4_processed, v9_processed)
rm(v4_processed, v9_processed)

# Convert to long format
myabundance <- myabundance %>%
  pivot_longer(
    cols = c("Araphid_pennate", "Raphid_pennate", "Coscinodiscophyceae", "Mediophyceae", "unknown_class"),
    names_to = "family",
    values_to = "sum_perc_euk"
  )


#### Plotting 


# Define the plotting function
plot_diatom_classes <- function(data, method) {
  # Filter data based on method and remove unknown_class
  tmp <- data[data$family != 'unknown_class' & data$method == method, ]
  
  # Set factor levels for family
  tmp$family <- factor(tmp$family, levels = c("Coscinodiscophyceae", "Mediophyceae", "Araphid_pennate", "Raphid_pennate"))
  
  # Create a presence column based on sum_perc_euk values
  tmp$presence <- ifelse(tmp$sum_perc_euk > 0, 'yes', 'no')
  
  # World map setup
  mapWorld <- borders("world", colour = NA, fill = "gray80")
  
  # Create map plot
  mp <- ggplot() + mapWorld + theme_bw() +
    geom_point(data = tmp[tmp$presence == 'yes', ], aes(y = event_latitude, x = event_longitude, size = sum_perc_euk), 
               shape = 1, color = 'blue', alpha = 0.5) +
    scale_size_area(name = '% euk\nreads', breaks = c(1, 5, 10, 20, 40, 80), max_size = 10) +
    facet_grid(family ~ method) +
    geom_point(data = tmp[tmp$presence == 'no', ], aes(y = event_latitude, x = event_longitude), 
               shape = 4, size = 2, color = 'red') +
    xlab('Longitude') + ylab('Latitude')
  
  # Set factor levels for size_fraction
  tmp$size_fraction <- factor(tmp$size_fraction, levels = c('180-2000', '20-180', '3/5-20', '0.8-5/2000'))
  
  # Create boxplot
  p <- ggplot(data = tmp, aes(x = size_fraction, y = sum_perc_euk)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_grid(family ~ method) +
    coord_flip() +
    scale_y_continuous(trans = 'log1p', breaks = c(1, 5, 10, 30, 60, 90)) +
    ylab('% eukaryotic reads') + xlab('Size fraction (μm)') +
    theme_bw()
  
  # Combine map and boxplot in a grid
  p_combined <- plot_grid(p, mp, ncol = 2, rel_widths = c(2, 3))
  
  # Display the combined plot
  print(p_combined)
  
  # Optionally save the plot to a file
#  ggsave(p_combined, file = paste0(method, '.pdf'), width = 8, height = 6)
}

# Run the plotting function
FigS7ab <- plot_diatom_classes(myabundance, 'V4')
FigS7ab

Fig4ab <- plot_diatom_classes(myabundance, 'V9')
Fig4ab

