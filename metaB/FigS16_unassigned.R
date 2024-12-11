
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
  
  # classify by assigned at genus level (genus/species) vs unassigned at genus level (only classified at phylum or class level)
  combined_data$genus <- gsub(x=combined_data$taxonomy, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Radial-centric-basal-Coscinodiscophyceae;', replacement = '')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Raphid-pennate;', replacement = '')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Polar-centric-Mediophyceae;', replacement = '')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Araphid-pennate;', replacement = '')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;Bacillariophyta_XX;Bacillariophyta_XXX;Bacillariophyta_XXX_sp.', replacement = 'unknown_genus')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = 'Root;Eukaryota;Stramenopiles;Ochrophyta;Bacillariophyta;Bacillariophyta_X;unclassified_Bacillariophyta_X', replacement = 'unknown_genus')
  combined_data$genus <- gsub(x=combined_data$genus, pattern = ';.*', replacement = '')
  combined_data[grep(x=combined_data$genus, pattern = 'X'),]$genus <- 'unknown_genus'
  combined_data[grep(x=combined_data$genus, pattern = 'unclassified'),]$genus <- 'unknown_genus'
  combined_data[combined_data$genus!='unknown_genus',]$genus <- 'known_genus'
  
  # ASVs per genus & class
  richness <- combined_data %>% count(genus, class, name = "row_count") 
  colnames(richness) <- c('genus', 'class', 'ASVs')

  # abundance sum per diatom genus
  myabundance <- combined_data %>% select(genus, class, samples$sample_id_pangaea)
  rm(combined_data)
  myabundance <- aggregate(data=myabundance, FUN=sum, . ~ genus + class)

  # reads per genus
  reads <- myabundance %>%
    select(genus, class, starts_with("TARA_")) %>%  # Select genus, class and TARA_ columns
    group_by(genus, class) %>%  # Group by genus and class
    summarise(across(starts_with("TARA_"), \(x) sum(x, na.rm = TRUE))) %>%  # Sum across TARA_ columns
    mutate(reads = rowSums(across(starts_with("TARA_"))))  # Calculate row sum per genus and class
  reads <- reads[, c("genus", "class", "reads")]

  richness_and_reads <- merge(x=richness, y=reads, all=T)
  
  # Convert to long format
  myabundance <- myabundance %>%
    pivot_longer(
      cols = starts_with("TARA_"),       # Select columns that start with "TARA_"
      names_to = "sample",               # Name for the new 'sample' column
      values_to = "counts"               # Name for the new 'counts' column
    )
  
  # merge genus counts with the total counts for normalization
  Indices <- merge(x=myabundance, y=Indices, by="sample", all.x=T)
  
  Indices$percEuk <- 100 * Indices$counts / Indices$ReadAbundance_Euks
  Indices$percDia <- 100 * Indices$counts / Indices$ReadAbundance_Diatoms
  
  # Merge with sample metadata
  Indices <- merge(Indices, samples, by.x="sample", by.y="sample_id_pangaea", all.x=TRUE)
  rm(samples)
  
  # calculate mean of technical replicates
  Indices <- aggregate(data=Indices, FUN=mean, cbind(percDia, percEuk) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction + genus + class)
  
  # Organize size fractions and mean equivalent size fractions
  Indices$size_fraction <- recode(Indices$size_fraction, "3-20" = "3/5-20", "5-20" = "3/5-20", "0.8-5" = "0.8-5/2000", "0.8->" = "0.8-5/2000", "0.8-20" = "0.8-5/2000")
  Indices$size_fraction <- factor(Indices$size_fraction, levels = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
  Indices <- aggregate(data=Indices, FUN=mean, cbind(percDia, percEuk) ~ 
                         event_latitude + event_longitude + ocean_region + station + depth + size_fraction + genus + class)
  
  #  Return both richness per genus and abundance per genus and sample:
  return(list(data_frame1 = Indices, data_frame2 = richness_and_reads))
}

# Process v4 and v9 data
v4_list <- process_diatom_data(v4, mysamples)
v4_processed <- v4_list$data_frame1
v4_processed$method <- 'V4'
v4_richness <- v4_list$data_frame2
rm(v4_list)

v9_list <- process_diatom_data(v9, mysamples)
v9_processed <- v9_list$data_frame1
v9_processed$method <- 'V9'
v9_richness <- v9_list$data_frame2
rm(v9_list)

myabundance <- rbind(v4_processed, v9_processed)
rm(v4_processed, v9_processed)
v4_richness$method <- 'V4'
v9_richness$method <- 'V9'
totals <- rbind(v4_richness, v9_richness)
rm(v4_richness, v9_richness)



#############################################

library(ggplot2)

totals_agg <- aggregate(data=totals, FUN=sum, cbind(ASVs, reads) ~ method + genus)

# Calculate the percentage of ASVs within each group
totals_agg <- totals_agg %>%
  group_by(method) %>%
  mutate(perc_ASVs = ASVs / sum(ASVs) * 100)

# Plot with percentage labels
p_ASVs <- ggplot(totals_agg, aes(y = ASVs, x = method, fill = genus)) + 
  geom_col() + 
  geom_text(aes(label = paste0(round(perc_ASVs, 1), "%")), 
            position = position_stack(vjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
) + 
  ggtitle("# ASVs") + 
  xlab("Method") + 
  ylab("# ASVs") +
  labs(fill = "Assignation")

p_ASVs


library(ggplot2)
# Calculate the percentage of reads within each group
totals_agg <- totals_agg %>%
  group_by(method) %>%
  mutate(perc_reads = reads / sum(reads) * 100)

# Plot with percentage labels
p_reads <- ggplot(totals_agg, aes(y = reads, x = method, fill = genus)) + 
  geom_col() + 
  geom_text(aes(label = paste0(round(perc_reads, 1), "%")), 
            position = position_stack(vjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
  legend.position = 'none') + 
  ggtitle("# reads") + 
  xlab("Method") + 
  ylab("# reads")

p_reads


library(cowplot)
p1 <- plot_grid(p_ASVs, p_reads, ncol = 2, rel_widths = c(1.5,1))

########

# Plot with percentage labels
p_ASVs_classes <- ggplot(totals, aes(y = ASVs, x = class, fill = genus)) + 
  geom_col(position = position_dodge(preserve = "single"), width = 0.6) +  # Dodge with equal width bars
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = 'none'
  ) + 
  ggtitle("# ASVs") + 
  xlab("taxon") + 
  ylab("# ASVs") +
  facet_grid(. ~ method) + 
  coord_flip()

p_ASVs_classes

# Plot with percentage labels
p_reads_classes <- ggplot(totals, aes(y = reads, x = class, fill = genus)) + 
  geom_col(position = position_dodge(preserve = "single"), width = 0.6) +  # Dodge with equal width bars
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = 'none'
  ) + 
  ggtitle("# reads") + 
  xlab("taxon") + 
  ylab("# reads") +
  facet_grid(. ~ method, scales = 'free') + 
  coord_flip()

p_reads_classes

p2 <- plot_grid(p_ASVs_classes, p_reads_classes, ncol = 2)


p00 <- plot_grid(p1, p2, ncol = 1, rel_heights = c(1,1.5))


########
myabundance$ocean_region <- gsub(x=myabundance$ocean_region, pattern = '].*', replacement = '')
myabundance$ocean_region <- gsub(x=myabundance$ocean_region, pattern = '\\[', replacement = '')

myabundance <- aggregate(data=myabundance, FUN=sum, cbind(percDia, percEuk) ~ station + depth + size_fraction + genus + ocean_region + method + event_latitude + event_longitude)


# Plot oceans
p_ocean <- ggplot(myabundance, aes(y = percEuk, x = ocean_region, fill = genus)) + 
  geom_boxplot(outlier.shape = 1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = 'none') + 
  scale_y_continuous(trans = scales::log1p_trans(),
                     breaks = c(1,10,20, 50,100),
                     labels = c(1,10,20, 50,100)) +
  facet_grid(genus ~ method) +
  ggtitle("Ocean") + 
  xlab("Ocean") + 
  ylab("% eukaryotic reads") +
  coord_flip() 

  p_ocean

mymap <- NULL
mapWorld <- borders("world", colour=NA, fill="gray90") # create a layer of borders
mymap <- ggplot() + mapWorld + theme_bw()
mymap <- mymap + geom_point(data=myabundance[myabundance$percEuk!=0,], aes(y=event_latitude, x=event_longitude, size=percEuk, color=genus), shape=16, alpha=0.5) + scale_size_area(max_size=10) + facet_grid(genus ~ method) 
mymap <- mymap + geom_point(data=myabundance[myabundance$percEuk==0,], aes(y=event_latitude, x=event_longitude), size=1, color="gray30", shape=4) + facet_grid(genus ~ method) + guides(size=guide_legend("% euk reads")) #+ labs(title = "map") 
mymap


myabundance$size_fraction <- factor(myabundance$size_fraction, levels=c('180-2000', '20-180', '3/5-20', '0.8-5/2000'))
  
  # Plot sizes
p_size <- ggplot(myabundance, aes(y = percEuk, x = size_fraction, fill = genus)) + 
    geom_boxplot(outlier.shape = 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          legend.position = 'none') + 
    scale_y_continuous(trans = scales::log1p_trans(),
                       breaks = c(1,10,20, 50,100),
                       labels = c(1,10,20, 50,100)) +
    facet_grid(genus ~ method) +
    ggtitle("Size fraction") + 
    xlab("Size fraction") + 
    ylab("% eukaryotic reads") +
    coord_flip() 
  
p_size

p01 <- plot_grid(p_size, p_ocean, ncol = 2, rel_widths = c(1,1.5))




#############################
library(reshape2)

unassigned <- myabundance[myabundance$genus=='unknown_genus',]
unassigned_wide <- dcast(unassigned, station + depth + size_fraction + ocean_region ~ method, value.var = "percEuk")
rm(unassigned)

# only samples with both V4 and V9
unassigned_wide <- unassigned_wide[!is.na(unassigned_wide$V4),]
unassigned_wide <- unassigned_wide[!is.na(unassigned_wide$V9),]

 unassigned_wide$size_fraction <- factor(unassigned_wide$size_fraction, levels=c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
 unassigned_wide$depth <- factor(unassigned_wide$depth, levels=c("SRF", "DCM"))

unassigned_wide$station <- gsub(unassigned_wide$station, pattern="^00", replacement="")
unassigned_wide$station <- gsub(unassigned_wide$station, pattern="^0", replacement="")
   

p <- ggplot(unassigned_wide, aes(y=V4, x=V9, color=ocean_region, label=station)) + 
		geom_text() +
		facet_grid(depth ~ size_fraction) +
		scale_y_continuous(trans = scales::log1p_trans(),
				breaks = c(1,10,20, 50,100),
				labels = c(1,10,20, 50,100), limits=c(0, 100)) +
		scale_x_continuous(trans = scales::log1p_trans(),
				breaks = c(1,10,20, 50,100),
				labels = c(1,10,20, 50,100), limits=c(0, 100)) +
		coord_fixed() +
		geom_abline(slope=1) + 
		 scale_color_manual(breaks = c("MS","RS","IO","SAO", "SO", "SPO","NPO","NAO","AO"),
		 values=c("darkblue","firebrick4","chocolate1","yellowgreen","black","lightblue3","yellow", "forestgreen","brown1")) +
    		xlab("unassigned V9 (% eukaryotic reads)") + ylab("unassigned V4 (% eukaryotic reads)") 
p
    
    
    
figureS16 <- plot_grid(p00, p01, p, ncol=1, rel_heights = c(1.5,1,1))
#ggsave(figureS16, file='FigS16.pdf')	

