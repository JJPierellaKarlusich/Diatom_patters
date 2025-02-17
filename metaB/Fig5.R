
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
  
  # ASVs per genus
    richness <- combined_data %>% count(genus, name = "row_count") 
    colnames(richness) <- c('genus', 'ASVs')
    richness <- richness[richness$genus!='unknown_genus',]
 
  # abundance sum per diatom genus
  myabundance <- combined_data %>% select(genus, class, samples$sample_id_pangaea)
  rm(combined_data)
  myabundance <- aggregate(data=myabundance, FUN=sum, . ~ genus + class)
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
  return(list(data_frame1 = Indices, data_frame2 = richness))
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
richness <- rbind(v4_richness, v9_richness)
rm(v4_richness, v9_richness)


##### genus #########################

genus_sum <- aggregate(data=myabundance, FUN=sum, cbind(percDia, percEuk) ~ genus + method)

genus_sum$genus <- factor(genus_sum$genus, levels=rev(unique(sort(genus_sum$genus))))

p_genus <- ggplot(genus_sum, aes(x=percEuk, y=genus, fill=method)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_x_continuous(trans = scales::log1p_trans(),
                     breaks = c(1,10,100,1000),
                     labels = c(1,10,100,1000)
  ) + 
  ggtitle("Abundance") + 
  xlab("% diatom reads") +  
  theme(legend.position = "none") +
  facet_grid(. ~ method) + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())

########### depth ####################

depth_sum <- aggregate(data=myabundance, FUN=sum, cbind(percDia, percEuk) ~ genus + method + depth)
depth_sum <- depth_sum[depth_sum$depth=='SRF' | depth_sum$depth=='DCM',]

depth_sum$genus <- factor(depth_sum$genus, levels=rev(unique(sort(depth_sum$genus))))

plot_depth <- ggplot(depth_sum, aes(x=percEuk, y=genus, fill=depth)) + 
    geom_bar(position="fill", stat="identity") + 
  ggtitle("Depth") + 
  xlab("% diatom reads") +  
  theme(legend.position = "none") +
  scale_fill_brewer(palette="Paired") + 
  scale_fill_manual(values = c("darkred", "lightgreen")) + 
  facet_grid(. ~ method) + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())


########### size ####################
size_sum <- aggregate(data=myabundance, FUN=sum, cbind(percDia, percEuk) ~ genus + method + size_fraction)
size_sum$size_fraction <- factor(size_sum$size_fraction, levels=c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))

size_sum$genus <- factor(size_sum$genus, levels=rev(unique(sort(size_sum$genus))))

plot_size <- ggplot(size_sum, aes(x=percEuk, y=genus, fill=size_fraction)) + 
     geom_bar(position="fill", stat="identity") + 
              ggtitle("Size") + theme(legend.position = "none") +
  scale_fill_manual(breaks = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"),
    values=c("#E31A1C","#1F78B4","#33A02C","purple")) +
   facet_grid(. ~ method) + 
  xlab("% diatom reads") +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())


########### occupancy ####################
occupancy <- myabundance[myabundance$percEuk!=0,]
occupancy <- unique(occupancy %>% select(genus, method, station))
occupancy <- aggregate(data=occupancy, FUN=length, . ~ genus + method)

occupancy$genus <- factor(occupancy$genus, levels=rev(unique(sort(occupancy$genus))))
  
plot_occupancy <- ggplot(occupancy, aes(x=station, y=genus, fill=method)) + 
    geom_bar(position="stack", stat="identity") + 
             ggtitle("Occupancy") + 
  xlab("# stations") +
  theme(legend.position = "none") +
   facet_grid(. ~ method) + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())

#####################

richness$genus <- factor(richness$genus, levels=rev(unique(sort(richness$genus))))

p_richness <- ggplot(richness, aes(x=ASVs, y=genus, fill=method)) + 
		geom_bar(position="dodge", stat="identity") + 
		scale_x_continuous(trans = scales::log1p_trans(),
				breaks = c(1,10,100,1000),
				labels = c(1,10,100,1000)
		) + 
		ggtitle("Richness") + 
  xlab("# ASVs") +
  theme(legend.position = "none") +
		facet_grid(. ~ method) + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())

######## classes ###############

classes <- unique(myabundance[c("genus", "class")])
classes$genus <- factor(classes$genus, levels=rev(unique(sort(classes$genus))))
classes$method <- ""

p_classes <- ggplot(classes, aes(x=1, y=genus, color=class, label=genus)) + 
  geom_text(hjust=1, size=3) +  # Align text to the left
  ggtitle("Genera") +  xlim(0.8, 1.005) +
  theme(legend.position = "none", 
        axis.text.y=element_blank(), 
        axis.title.y=element_blank(),
        axis.title.x = element_text(color = "white"),
        axis.text.x = element_text(color = "white"),  # Set x-axis labels to white
        panel.background = element_rect(fill = "white", color = NA),  # White panel background
        plot.background = element_rect(fill = "white", color = NA),   # White plot background
        panel.grid.major = element_line(color = "gray90"),  # Light grid lines
        panel.grid.minor = element_blank()) +  # No minor grid lines
  facet_grid(. ~ method)  +
  scale_color_manual(values = c("Mediophyceae" = "DarkBlue",
                                "Coscinodiscophyceae" = "LightBlue",
                                "Raphid_pennate" = "DarkGreen",
                                "Araphid_pennate" = "Orange")) 


#######################

library(cowplot)
p <- plot_grid(p_classes, p_richness, p_genus, plot_size, plot_depth, plot_occupancy, ncol=6, rel_widths=c(1,1.1,1.1,1.1,1.1,1.1 ))
p
#ggsave(p, file="figure_summary_May2024.pdf", height=12, width=10)
