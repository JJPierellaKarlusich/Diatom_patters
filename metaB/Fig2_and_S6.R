
library(data.table)
library(tidyverse)
library(ggplot2)
library(vegan)
library(cowplot)


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
  combined_data <- combined_data %>%
    select(amplicon, taxonomy, confidence, total, spread, sequence, samples$sample_id_pangaea)
  
  combined_data$total_subset <- rowSums(combined_data %>% select(samples$sample_id_pangaea))

  # ASVs with at least 3 counts in the whole selected dataset
  combined_data$total_subset <- rowSums(combined_data %>% select(samples$sample_id_pangaea))
  combined_data <- combined_data[combined_data$total_subset >= 3,]
  combined_data$total_subset <- NULL
  
  # ASVs with ocurrance in at least 2 different samples the whole selected 
  sample_columns <- combined_data[, 7:ncol(combined_data)] #Identify the sample columns (7th to the last)
  occurrences <- rowSums(sample_columns > 0) #Calculate the number of samples where each species occurs
  rm(sample_columns)
  combined_data <- combined_data[occurrences >= 2, ] #Filter rows with occurrences in at least two different samples
  rm(occurrences)
    
  myabundance <- combined_data %>% select(samples$sample_id_pangaea)

  # Compute diatom read abundance
  myabundance <- t(myabundance)
  ReadAbundance <- rowSums(myabundance)
  ReadAbundance <- as.data.frame(ReadAbundance)
  colnames(ReadAbundance) <- "readsDiatoms"
  
  # Calculate diatom diversity metrics
  Shannon = diversity(myabundance, index="shannon",MARGIN = 1,base = exp(1))
  Shannon <- as.data.frame(Shannon)
  colnames(Shannon) <- "Shannon"
  
  expShannon <- exp(Shannon)
  expShannon <- as.data.frame(expShannon)
  colnames(expShannon) <- "expShannon"
  
  Indices <- merge(x=Shannon, y=expShannon, by="row.names", all=T)
  colnames(Indices)[1] <- "sample"
  rm(Shannon, expShannon)
  Indices <- merge(x=Indices, y=ReadAbundance, by.x='sample', by.y='row.names')

  # Calculate total eukaryotic read abundance and normalize diatom relative abundance
  euk_data <- data %>% select(samples$sample_id_pangaea)
  euk_data <- t(euk_data)
  ReadAbundance <- rowSums(euk_data, na.rm = TRUE)
  rm(euk_data)
  ReadAbundance <- data.frame(ReadAbundance_Euks = ReadAbundance)
  Indices <- merge(x=Indices, y=ReadAbundance, by.x="sample", by.y="row.names", all=T)
  
  Indices$percEuk <- 100 * Indices$readsDiatoms / Indices$ReadAbundance_Euks
  
  # Merge with sample metadata
  Indices <- merge(Indices, samples, by.x="sample", by.y="sample_id_pangaea", all.x=TRUE)

  # calculate mean of technical replicates
  Indices <- aggregate(data=Indices, FUN=mean, cbind(percEuk, readsDiatoms, Shannon, expShannon) ~ 
    event_latitude + event_longitude + ocean_region + station + depth + size_fraction + abs_lat)
  
  # Organize size fractions and mean equivalent size fractions
  Indices$size_fraction <- recode(Indices$size_fraction, "3-20" = "3/5-20", "5-20" = "3/5-20", "0.8-5" = "0.8-5/2000", "0.8->" = "0.8-5/2000", "0.8-20" = "0.8-5/2000")
  Indices$size_fraction <- factor(Indices$size_fraction, levels = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
  Indices <- aggregate(data=Indices, FUN=mean, cbind(percEuk, readsDiatoms, Shannon, expShannon) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction + abs_lat)
  
  
}



# Process v4 and v9 data
v4_processed <- process_diatom_data(v4, mysamples)
v4_processed$method <- 'V4'
v9_processed <- process_diatom_data(v9, mysamples)
v9_processed$method <- 'V9'
rm(process_diatom_data)

#############################
#############################

# Define function for plotting diatom latitudinal gradients and biogeography
ploteador <- function(indices){

p_abundance_lat <- ggplot(data=indices, aes(x=event_latitude, y=percEuk)) + 
  geom_point(shape=1) + geom_smooth() +
  facet_grid(size_fraction ~ method) + 
  coord_flip() + 
  ylab('Diatom relative abundance\n(% eukaryotic reads)') + xlab('Latitude')

p_ashannon_lat <- ggplot(data=indices, aes(x=event_latitude, y=expShannon)) + 
  geom_point(shape=1) + geom_smooth() +
  facet_grid(size_fraction ~ method) + 
  coord_flip() + 
  ylab('Diatom exp Shannon') + xlab('Latitude')

library("wesanderson")
pal <- wes_palette("Zissou1", 21, type = "continuous")

tmp <- indices[indices$depth=="SRF",]
mymap <- NULL
mapWorld <- borders("world", colour=NA, fill="gray80") # create a layer of borders
mymap <- ggplot() + mapWorld + theme_bw()
mymap <- mymap + geom_point(data=tmp[tmp$percEuk!=0,], aes(y=event_latitude, x=event_longitude, size=percEuk, color=expShannon), shape=16) + scale_size_area(max_size=10) + facet_grid(size_fraction ~ method) +  scale_color_gradientn(colours=pal)
mymap <- mymap + geom_point(data=tmp[tmp$percEuk==0,], aes(y=event_latitude, x=event_longitude), size=1, color="gray30", shape=4) + facet_grid(size_fraction ~ method) + guides(size=guide_legend("% euk reads")) #+ labs(title = "map") 

plot_grid(mymap, p_abundance_lat, p_ashannon_lat, ncol=3, rel_widths = c(2,1,1))

}

########

v4_plot <- ploteador(v4_processed)
v4_plot
v9_plot <- ploteador(v9_processed)
v9_plot

