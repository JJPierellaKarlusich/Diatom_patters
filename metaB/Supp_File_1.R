##############

# Load required library

library(data.table)
library(tidyverse)
library(ggplot2)
library(vegan)
library(cowplot)
library(tidyr)
library(dplyr)


#### Formating the data

# Load v4 and v9 tables
# ASVs with taxonomy with confidence > 50% by idtaxa

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
  
  # Select only relevant samples 
  samples <- samples[samples$sample_id_pangaea %in% colnames(combined_data),]
  combined_data <- combined_data %>% select(amplicon, taxonomy, confidence, total, spread, sequence, samples$sample_id_pangaea)
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
  
  # Calculate total diatom read abundance for normalizing diatom genera relative abundance (and compare with the normalization by total eukaryotes)
  dia_data <- combined_data %>% select(samples$sample_id_pangaea)
  dia_data <- t(dia_data)
  dia_data <- rowSums(dia_data, na.rm = TRUE)
  dia_data <- data.frame(ReadAbundance_Diatoms = dia_data)
  
  # Calculate total eukaryotic read abundance for normalizing diatom genera relative abundance
  euk_data <- data %>% select(samples$sample_id_pangaea)
  euk_data <- t(euk_data)
  euk_data <- rowSums(euk_data, na.rm = TRUE)
  euk_data <- data.frame(ReadAbundance_Euks = euk_data)
  
  # merge total counts of diatoms and eukaryotes
  Indices <- merge(x=dia_data, y=euk_data, by="row.names", all.x=T)
  rm(dia_data, euk_data)
  colnames(Indices)[1] <- 'sample'
  
  # classify by diatom class
  #"Mediophyceae" # polar and 'non-polar' centric
  #"Coscinodiscophyceae" #Radial-centric-basal
  combined_data$class <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Radial-centric-basal-Coscinodiscophyceae.*', replacement = 'Coscinodiscophyceae')
  combined_data$class <- gsub(x=combined_data$class, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Raphid-pennate.*', replacement = 'Raphid_pennate')
  combined_data$class <- gsub(x=combined_data$class, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Polar-centric-Mediophyceae.*', replacement = 'Mediophyceae')
  combined_data$class <- gsub(x=combined_data$class, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Araphid-pennate.*', replacement = 'Araphid_pennate')
  combined_data$class <- gsub(x=combined_data$class, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Bacillariophyta_XX;Bacillariophyta_XXX;Bacillariophyta_XXX_sp.', replacement = 'unknown_class')
  combined_data$class <- gsub(x=combined_data$class, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;unclassified_Bacillariophyta_X', replacement = 'unknown_class')
  
  # diatom genera
  combined_data$genus <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Radial-centric-basal-Coscinodiscophyceae;', replacement = '')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Raphid-pennate;', replacement = '')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Polar-centric-Mediophyceae;', replacement = '')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Araphid-pennate;', replacement = '')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Bacillariophyta_XX;Bacillariophyta_XXX;Bacillariophyta_XXX_sp.', replacement = 'unknown_genus')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;unclassified_Bacillariophyta_X', replacement = 'unknown_genus')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = ';.*', replacement = '')
  combined_data[grep(x=combined_data$genus, pattern = 'X'),]$genus <- 'unknown_genus'
  combined_data[grep(x=combined_data$genus, pattern = 'unclassified'),]$genus <- 'unknown_genus'
  
  combined_data[combined_data$genus=='Odontella' | combined_data$genus=='Trieres',]$genus <- 'Trieres/Odontella'
  
  # remove unassigned ASVs
  combined_data <- combined_data[combined_data$genus!='unknown_genus',]
  
  
  # merge class and genus
  combined_data$genus <- paste(combined_data$genus, combined_data$class, sep=' | ')

  # matrix of ASVs (named by genus) vs samples
  combined_data <- combined_data %>% select(genus, samples$sample_id_pangaea)

  # read abundance sum per diatom genus
  myabundance <- aggregate(data=combined_data, FUN=sum, . ~ genus)
  myabundance <- myabundance[myabundance$genus!='unknown_genus',]
  
  # Convert to long format
  myabundance <- myabundance %>%
    pivot_longer(
      cols = starts_with("TARA_"),       # Select columns that start with "TARA_"
      names_to = "sample",               # Name for the new 'sample' column
      values_to = "counts"               # Name for the new 'counts' column
    )
  
  # merge genus counts with the total counts for normalization
  Indices <- merge(x=myabundance, y=Indices, by="sample", all.x=T)
  rm(myabundance)
  
  Indices$percEuk <- 100 * Indices$counts / Indices$ReadAbundance_Euks
  Indices$percDia <- 100 * Indices$counts / Indices$ReadAbundance_Diatoms
  
  # Merge with sample metadata
  Indices <- merge(Indices, samples, by.x="sample", by.y="sample_id_pangaea", all.x=TRUE)
  rm(samples)
  
  # calculate mean of technical replicates
  Indices <- aggregate(data=Indices, FUN=mean, cbind(percDia, percEuk) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction + genus)
  
  # Organize size fractions and mean equivalent size fractions
  Indices$size_fraction <- recode(Indices$size_fraction, "3-20" = "3/5-20", "5-20" = "3/5-20", "0.8-5" = "0.8-5/2000", "0.8->" = "0.8-5/2000", "0.8-20" = "0.8-5/2000")
  Indices$size_fraction <- factor(Indices$size_fraction, levels = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
  Indices <- aggregate(data=Indices, FUN=mean, cbind(percDia, percEuk) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction + genus)
  
  # Return both Indices (noromalized read abundance per genus and sample) and myabundance (community matrix) as a list
  return(list(Indices = Indices, combined_data = combined_data))
}


# Process v4 data
v4_processed_list <- process_diatom_data(v4, mysamples)
v4_processed <- v4_processed_list$Indices
v4_processed$method <- 'V4'
v4_matrix <- v4_processed_list$combined_data
rm(v4_processed_list)

# Process v9 data
v9_processed_list <- process_diatom_data(v9, mysamples)
v9_processed <- v9_processed_list$Indices
v9_processed$method <- 'V9'
v9_matrix <- v9_processed_list$combined_data
rm(v9_processed_list)

myabundance <- rbind(v4_processed, v9_processed)
rm(v4_processed, v9_processed)
rm(process_diatom_data)

######### calculate shannon index for each genus #################
library(vegan)

#function
shannon_calculator <- function(GENUS, MATRIX){
  barcodes_genus <- MATRIX
  barcodes_genus <- barcodes_genus[barcodes_genus$genus==GENUS,]
  barcodes_genus$genus <- NULL
  barcodes_genus <- t(barcodes_genus)
  
  ReadAbundance <- rowSums(barcodes_genus)
  ReadAbundance <- as.data.frame(ReadAbundance)
  colnames(ReadAbundance) <- "ReadAbundance"
  
  ASVnumber <- specnumber(barcodes_genus)
  ASVnumber <- as.data.frame(ASVnumber)
  colnames(ASVnumber) <- "ASVnumber"
  
  tmp1 <- merge(x=ReadAbundance, y=ASVnumber, by="row.names", all=T)
  
  Shannon = diversity(barcodes_genus,index="shannon",MARGIN = 1,base = exp(1))
  Shannon <- as.data.frame(Shannon)
  colnames(Shannon) <- "Shannon"
  
  expShannon <- exp(Shannon)
  expShannon <- as.data.frame(expShannon)
  colnames(expShannon) <- "expShannon"
  
  tmp2 <- merge(x=Shannon, y=expShannon, by="row.names", all=T)
  
  Indices <- merge(x=tmp1, y=tmp2, by="Row.names", all=T)
  colnames(Indices)[1] <- "sample"
  
  #Indices = data.frame(ReadAbundance,ASVnumber,Shannon,expShannon)
  
  rm(tmp1, tmp2, barcodes_genus, ReadAbundance,ASVnumber,Shannon,expShannon)
  
  Indices$genus <- GENUS
  return(Indices)
}

######## v4 shannon ######## 
genus_ids_v4 <- unique(v4_matrix$genus) # Get unique genus IDs
genus_indices_v4 <- data.frame() # Initialize an empty data frame to store results
# Loop through each genus and calculate Shannon indices
for (genus in genus_ids_v4) {
  genus_index_v4 <- shannon_calculator(genus, v4_matrix)   # Calculate the Shannon index for the current genus
  genus_indices_v4 <- rbind(genus_indices_v4, genus_index_v4)  # Append the result to the genus_indices data frame
}
rm(genus_ids_v4, v4_matrix)
genus_indices_v4$method <- 'V4'

######## v9 shannon ######## 
genus_ids_v9 <- unique(v9_matrix$genus) # Get unique genus IDs
genus_indices_v9 <- data.frame() # Initialize an empty data frame to store results
# Loop through each genus and calculate Shannon indices
for (genus in genus_ids_v9) {
  genus_index_v9 <- shannon_calculator(genus, v9_matrix)   # Calculate the Shannon index for the current genus
  genus_indices_v9 <- rbind(genus_indices_v9, genus_index_v9)  # Append the result to the genus_indices data frame
}
rm(genus_ids_v9, v9_matrix)
genus_indices_v9$method <- 'V9'

shannon <- rbind(genus_indices_v4, genus_indices_v9)
rm(genus_indices_v4, genus_indices_v9)
rm(shannon_calculator)

# Merge with sample metadata
shannon <- merge(x=shannon, y=mysamples, by.x="sample", by.y="sample_id_pangaea", all.x=TRUE)
#rm(mysamples)

# calculate mean of technical replicates
shannon <- aggregate(data=shannon, FUN=mean, cbind(ReadAbundance, ASVnumber, Shannon, expShannon) ~ 
                       event_latitude + event_longitude + ocean_region + station + depth + size_fraction + genus + method)

# Organize size fractions and mean equivalent size fractions
shannon$size_fraction <- recode(shannon$size_fraction, "3-20" = "3/5-20", "5-20" = "3/5-20", "0.8-5" = "0.8-5/2000", "0.8->" = "0.8-5/2000", "0.8-20" = "0.8-5/2000")
shannon$size_fraction <- factor(shannon$size_fraction, levels = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
shannon <- aggregate(data=shannon, FUN=mean, cbind(ReadAbundance, ASVnumber, Shannon, expShannon) ~ 
                       event_latitude + event_longitude + ocean_region + station + depth + size_fraction + genus + method)

### merge Shannon and normalized read abundance ####
myabundance <- merge(x=myabundance, y=shannon, all=T, by=c("event_latitude", "event_longitude", "ocean_region", "station", "depth", "size_fraction", "genus", "method"))
rm(shannon)

###################################################################
######## Plot biogeography ########################## 
###################################################################

#####################################

library("wesanderson")
pal <- wes_palette("Zissou1", 21, type = "continuous")

mapeador <- function(MYGENUS) {
  tmp <- myabundance[myabundance$genus==MYGENUS & myabundance$depth=="SRF",]
  mp1 <- NULL
  mapWorld <- borders("world", colour=NA, fill="gray80") # create a layer of borders
  mp1 <- ggplot() + mapWorld + theme_bw()
  mp1 <- mp1 + geom_point(data=tmp[tmp$percEuk!=0,], aes(y=event_latitude, x=event_longitude, size=percEuk, color=expShannon), shape=16) + scale_size_area(max_size=10) + facet_grid(size_fraction ~ method) +  scale_color_gradientn(colours=pal)
  mp1 <- mp1 + geom_point(data=tmp[tmp$percEuk==0,], aes(y=event_latitude, x=event_longitude), size=1, color="gray30", shape=4) + facet_grid(size_fraction ~ method) + labs(title = MYGENUS) + guides(size=guide_legend("% euk reads"))
  return(mp1)
  rm(tmp, mp1)
}


# Get unique genus IDs
genus_list <- unique(myabundance$genus)

# Initialize an empty list to store results
mapped_results <- list()

# Loop through each genus and apply the mapeador function
for (genus in genus_list) {
  print(genus)
  # Construct the argument for mapeador
  # Apply the mapeador function
  mapped_results[[genus]] <- mapeador(genus)
  genusID <- gsub(genus, pattern = ' .*', replacement = '')
  print(genusID)
  ggsave(mapped_results[[genus]], file=paste(genusID, '.png', sep=''))
  
}

# View the results
mapped_results

# Chaetoceros <- mapeador('Chaetoceros | Mediophyceae')

