library(tidyverse)
library(ggplot2)
library(cowplot)
library(plsdepot)
library(ggplotify)
library(tidyr)
library(reshape2)
library(Hmisc)
library(GGally)


################################################################################
###### Calculate diatom relative abundances ######
################################################################################

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
  
  # Select only relevant samples 
  samples <- samples[samples$sample_id_pangaea %in% colnames(combined_data),]
  combined_data <- combined_data %>%
    select(amplicon, taxonomy, confidence, total, spread, sequence, samples$sample_id_pangaea)
  
  combined_data$total_subset <- rowSums(combined_data %>% select(samples$sample_id_pangaea))
  combined_data <- combined_data[combined_data$total_subset != 0,]
  
  myabundance <- combined_data %>% select(samples$sample_id_pangaea)
  
  # Compute diatom read abundance
  myabundance <- t(myabundance)
  ReadAbundance_diatoms <- rowSums(myabundance)
  ReadAbundance_diatoms <- as.data.frame(ReadAbundance_diatoms)
  colnames(ReadAbundance_diatoms) <- "readsDiatoms"
  
  # Calculate total eukaryotic read abundance and normalize diatom relative abundance
  euk_data <- data %>% select(samples$sample_id_pangaea)
  euk_data <- t(euk_data)
  ReadAbundance <- rowSums(euk_data, na.rm = TRUE)
  rm(euk_data)
  ReadAbundance <- data.frame(ReadAbundance_Euks = ReadAbundance)
  Indices <- merge(x=ReadAbundance_diatoms, y=ReadAbundance, by="row.names", all=T)
  
  Indices$percEuk <- 100 * Indices$readsDiatoms / Indices$ReadAbundance_Euks
  
  # Merge with sample metadata
  Indices$sample <- Indices$Row.names
  Indices$Row.names <- NULL
  Indices <- merge(Indices, samples, by.x="sample", by.y="sample_id_pangaea", all.x=TRUE)
  
  # calculate mean of technical replicates
  Indices <- aggregate(data=Indices, FUN=mean, cbind(percEuk, readsDiatoms) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction + biomeplot)
  
  # Organize size fractions and mean equivalent size fractions
  Indices$size_fraction <- recode(Indices$size_fraction, "3-20" = "3/5-20", "5-20" = "3/5-20", "0.8-5" = "0.8-5/2000", "0.8->" = "0.8-5/2000", "0.8-20" = "0.8-5/2000")
  Indices$size_fraction <- factor(Indices$size_fraction, levels = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
  Indices <- aggregate(data=Indices, FUN=mean, cbind(percEuk, readsDiatoms) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction + biomeplot)
  
  
}



# Process v4 and v9 data
v4_processed <- process_diatom_data(v4, mysamples)
v4_processed$method <- 'V4'
v9_processed <- process_diatom_data(v9, mysamples)
v9_processed$method <- 'V9'

v4v9 <- rbind(v4_processed, v9_processed)
rm(v4_processed, v9_processed)

v4v9$station <- gsub(x=v4v9$station, pattern = '^00', replacement = '')
v4v9$station <- gsub(x=v4v9$station, pattern = '^0', replacement = '')

rm(mysamples, v4, v9, process_diatom_data)

#########################################################################################################
########### plot ################
#########################################################################################################


# Ocean

v4v9$ocean <- gsub(v4v9$ocean_region, pattern="\\[", replacement="")
v4v9$ocean <- gsub(v4v9$ocean, pattern="\\].*", replacement="")
v4v9$ocean <- factor(v4v9$ocean, levels=c("NAO","MS", "RS", "IO","SAO","SO","SPO","NPO","AO"))

p_OR <- ggplot(data=v4v9, aes(x=ocean, y=percEuk)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=1, notch=FALSE) + 
  facet_grid(. ~ method) + 
  coord_flip() + 
  ylab("Diatom relative abundance (% eukaryotic reads)") + 
  EnvStats::stat_n_text(y.pos = 94, size=3)
p_OR


# water later

p_depth <- ggplot(data=v4v9, aes(x=depth, y=percEuk)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=1.5, notch=FALSE) + 
  facet_grid(. ~ method) + 
  coord_flip() + EnvStats::stat_n_text(y.pos = 112, size=3) + 
  ylab("Diatom relative abundance (% eukaryotic reads)") + 
  ylim(0,120)

p_depth

# biome

p_biome <- ggplot(data=v4v9, aes(x=biomeplot, y=percEuk)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=1.5, notch=FALSE) + 
  facet_grid(. ~ method) + 
  coord_flip() + EnvStats::stat_n_text(y.pos = 112, size=3) + 
  ylab("Diatom relative abundance (% eukaryotic reads)") + 
  ylim(0,120)

p_biome

# size fraction

v4v9$size_fraction <- factor(v4v9$size_fraction, levels=c("180-2000", "20-180", "3/5-20", "0.8-5/2000"))

p_size <- ggplot(data=v4v9, aes(x=size_fraction, y=percEuk)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=1, outlier.size=1, notch=FALSE) + 
  facet_grid(. ~ method) + coord_flip() + 
  ylab("Diatom relative abundance (% eukaryotic reads)") + 
  EnvStats::stat_n_text(y.pos = 112, size=3) + 
  ylim(0, 120)

p_size

FigS5 <- plot_grid(p_size, p_depth, p_OR, p_biome, ncol=1) #, ncol=3, rel_widths=c(2,1,1))
FigS5

