# Load required library

library(data.table)
library(tidyverse)
library(ggplot2)
library(vegan)
library(cowplot)
library(tidyr)
library(dplyr)


####################################################### metaB #####################################
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
  
  # abundance sum per diatom genus
  myabundance <- combined_data %>% select(genus, class, samples$sample_id_pangaea)
  rm(combined_data)
  myabundance <- aggregate(data=myabundance, FUN=sum, . ~ genus + class)
  
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
  
  #  Return abundance per genus and sample:
  return(Indices)
}



# Process v4 and v9 data
v4_processed <- process_diatom_data(v4, mysamples)
v4_processed$method <- 'V4'

v9_processed <- process_diatom_data(v9, mysamples)
v9_processed$method <- 'V9'

myabundance <- rbind(v4_processed, v9_processed)
rm(v4_processed, v9_processed)
myabundance <- myabundance[myabundance$genus!='unknown_genus',]


metaB_microsize <-  myabundance[myabundance$size_fraction=='20-180',]
rm(myabundance, v4, v9, mysamples, process_diatom_data)


# Due to similar morphologies, we merged the counts for Actinocyclus & Coscinodiscus (and centric diatoms >30 µm in length).
metaB_microsize[
  metaB_microsize$genus=="Actinocyclus" |  
    metaB_microsize$genus=="Coscinodiscus" ,]$genus <- "Actinocyclus/Coscinodiscus/centric>30µm"

metaB_microsize <- aggregate(data=metaB_microsize, FUN=sum, 
                             cbind(percDia, percEuk) ~ 
                               station + depth + size_fraction + event_latitude + event_longitude + ocean_region + genus + class + method)

classes <- unique(metaB_microsize[c("genus", "class")])


##################################################### microscopy ##################################

# upload microscopy data
microscopy_microsize <- read.csv(file="datasets/Microscopy/diatom_microscopy.tsv", sep="\t", header=T)

######## standarize genera for comparisson with metaB:
microscopy_microsize[microscopy_microsize$genus=="Shionodiscus",]$genus <- "Thalassiosira"
microscopy_microsize[microscopy_microsize$genus=="Odontella",]$genus <- "Trieres/Odontella"
# Due to similar morphologies, we merged the counts for Actinocyclus, Coscinodiscus and centric diatoms >30 µm in length
microscopy_microsize[microscopy_microsize$genus=="Actinocyclus" | microscopy_microsize$genus=="Coscinodiscus" | microscopy_microsize$genus=="Centric diatoms 30-50 µm " | microscopy_microsize$genus=="Centric 60 µm",]$genus <- "Actinocyclus/Coscinodiscus/centric>30µm"
# We also merged Fossulaphycus (previously Fossula; never cultivated or sequenced) with Fragilaria due to 
# similar morphologies and co occurrence (microscopy detection of Fossulaphycus only at one Arctic station, 
# where V4 detected a high abundance of Fragilaria reads,
microscopy_microsize[microscopy_microsize$genus=='Fossula',]$genus <- 'Fragilaria'
# remove unassigned
microscopy_microsize <- microscopy_microsize[microscopy_microsize$genus!="Centric 10 µm",]
microscopy_microsize <- microscopy_microsize[microscopy_microsize$genus!="Pennate",]
# sum equivalent genera (eg., Actinocyclus + Coscinodiscus + Centric diatoms 30-50 µm + Centric diatoms 60 µm )
microscopy_microsize <- aggregate(data=microscopy_microsize, FUN=sum, percDia ~ station + depth + genus + method)

################ merge metaB & microscopy

# standarize both dataframes
microscopy_microsize <- microscopy_microsize %>% select('station', 'depth', 'genus', 'percDia', 'method')
metaB_microsize <- metaB_microsize %>% select('station', 'depth', 'genus', 'percDia', 'method')
metaB_microsize$station <- gsub(metaB_microsize$station, pattern = '^00', replacement = '')
metaB_microsize$station <- gsub(metaB_microsize$station, pattern = '^0', replacement = '')

# Get the unique station-depth pairs from microscopy_microsize
unique_st_depth <- unique(microscopy_microsize[c("station", "depth")])
# Filter metaB_microsize based on these pairs
metaB_microsize <- merge(x=unique_st_depth, y=metaB_microsize, by=c("station", "depth"), all=F)
rm(unique_st_depth)

# Get the unique station-depth pairs from metaB_microsize
v4 <- metaB_microsize[metaB_microsize$method=='V4',]
v9 <- metaB_microsize[metaB_microsize$method=='V9',]
unique_st_depth_v4 <- unique(v4[c("station", "depth")])
unique_st_depth_v9 <- unique(v9[c("station", "depth")])
rm(v4,v9)
# Filter metaB_microsize based on these pairs
microscopy_microsize <- merge(x=unique_st_depth_v4, y=microscopy_microsize, by=c("station", "depth"), all.x=F)
microscopy_microsize <- merge(x=unique_st_depth_v9, y=microscopy_microsize, by=c("station", "depth"), all.x=F)
rm(unique_st_depth_v4, unique_st_depth_v9)

# remove samples present in microscopy but not in metaB:
#station depth genus percDia method
#079   DCM  <NA>      NA   <NA>
#125   MIX  <NA>      NA   <NA>

# merge metaB & microscopy
mydata <- rbind(metaB_microsize, microscopy_microsize)
rm(metaB_microsize, microscopy_microsize)


#################### top 20 ########################
# sum counts per genus
genus_sum <- aggregate(data=mydata, FUN=sum, percDia ~ genus + method)

genus_sum[genus_sum$percDia==0,]
genus_sum <- genus_sum[genus_sum$percDia!=0,]

MicroscopyrankDia <- genus_sum[genus_sum$method=='microscopy',]
MicroscopyrankDia$normalization <- 'diatoms'
MicroscopyrankDia$rank <- rank(-MicroscopyrankDia$percDia)

V4rankDia <- genus_sum[genus_sum$method=='V4',]
V4rankDia$normalization <- 'diatoms'
V4rankDia$rank <- rank(-V4rankDia$percDia)

V9rankDia <- genus_sum[genus_sum$method=='V9',]
V9rankDia$normalization <- 'diatoms'
V9rankDia$rank <- rank(-V9rankDia$percDia)

genus_sum <- rbind(MicroscopyrankDia, V4rankDia, V9rankDia)
rm(MicroscopyrankDia, V4rankDia, V9rankDia)


top20 <- genus_sum[genus_sum$rank<=20,]



######################### #########################
######################### #########################
######################### #########################

##### Top 20 genus #########################

genus_sum_top20 <- genus_sum[genus_sum$genus %in%  unique(top20$genus),]


# ranking for plotting in the same order
rankings <- genus_sum_top20
rankings$rank_for_plotting <-  paste(rankings$method, rankings$normalization, sep='_')
# reshape the data into wide format
rankings_wide <- rankings %>%
  select(genus, rank_for_plotting, rank) %>%
  pivot_wider(
    names_from = rank_for_plotting,
    values_from = rank
  )
rankings_wide <- as.data.frame(rankings_wide)
rankings_wide[is.na(rankings_wide$V4_diatoms),]$V4_diatoms <- 21
rankings_wide[is.na(rankings_wide$V9_diatoms),]$V9_diatoms <- 21
rankings_wide[is.na(rankings_wide$microscopy_diatoms),]$microscopy_diatoms <- 21

rankings_wide[rankings_wide$V4_diatoms>20,]$V4_diatoms <- 21
rankings_wide[rankings_wide$V9_diatoms>20,]$V9_diatoms <- 21
rankings_wide[rankings_wide$microscopy_diatoms>20,]$microscopy_diatoms <- 21

# Arrange the dataframe by the priority columns in descending order
rankings_wide <- rankings_wide %>%
  arrange(
    microscopy_diatoms,
    V4_diatoms, 
    V9_diatoms
  ) %>%
  mutate(order_rank = row_number())  # Add a column with ranks from 1 to 29 based on the ordering



# order genera by ranking 
genus_sum_top20$genus <- factor(genus_sum_top20$genus, levels=rankings_wide[order(-rankings_wide$order_rank),]$genus)





##### genus #########################


p_genus_dia <- ggplot(genus_sum_top20, aes(x=percDia, y=genus, fill=method)) + 
  geom_bar(position="dodge", stat="identity") + 
  scale_x_continuous(trans = scales::log1p_trans(),
                     breaks = c(1,10,100,1000),
                     labels = c(1,10,100,1000)
  ) + 
  ggtitle("Abundance normalized\nby diatoms") + 
  xlab("sum normalized reads") +  
  theme(legend.position = "none") +
  facet_grid(. ~ method) + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())


########### occupancy ####################
occupancy <- mydata[mydata$percDia!=0,]
occupancy <- unique(occupancy %>% select(genus, method, station))
occupancy <- aggregate(data=occupancy, FUN=length, . ~ genus + method)

occupancy_top20 <- occupancy[occupancy$genus %in%  unique(top20$genus),]
rm(occupancy)

# order genera by ranking in V4
occupancy_top20$genus <- factor(occupancy_top20$genus, levels=rankings_wide[order(-rankings_wide$order_rank),]$genus)

plot_occupancy <- ggplot(occupancy_top20, aes(x=station, y=genus, fill=method)) + 
  geom_bar(position="stack", stat="identity") + 
  ggtitle("Occupancy\n ") + 
  xlab("# stations") +
  theme(legend.position = "none") +
  facet_grid(. ~ method) + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank())


######## classes ###############

#classes <- unique(metaB_microsize[c("genus", "class")])
classes_top20 <- classes[classes$genus %in% unique(top20$genus),]
rm(classes)

# order genera by ranking in V4
classes_top20$genus <- factor(classes_top20$genus, levels=rankings_wide[order(-rankings_wide$order_rank),]$genus)
classes_top20$method <- ""

p_classes <- ggplot(classes_top20, aes(x=1, y=genus, color=class, label=genus)) + 
  geom_text(hjust=1, size=3) +  # Align text to the left
  ggtitle("Genera\n ") +  xlim(0.8, 1.005) +
  theme(legend.position = "none", 
        axis.text.y=element_blank(), 
        axis.title.y=element_blank(),
        axis.title.x = element_text(color = "white"),
        axis.text.x = element_text(color = "white"),  # Set x-axis labels to white
        panel.background = element_rect(fill = "white", color = NA),  # White panel background
        plot.background = element_rect(fill = "white", color = NA),   # White plot background
        panel.grid.major = element_line(color = "white"),  # Light grid lines
        panel.grid.minor = element_blank()) +  # No minor grid lines
  facet_grid(. ~ method)  +
  scale_color_manual(values = c("Mediophyceae" = "DarkBlue",
                                "Coscinodiscophyceae" = "LightBlue",
                                "Raphid_pennate" = "DarkGreen",
                                "Araphid_pennate" = "Orange")) 

p_classes
########### Rank ####################


# order genera by ranking in V4
genus_sum_top20$genus <- factor(genus_sum_top20$genus, levels=rankings_wide[order(-rankings_wide$order_rank),]$genus)

# different color if not top20 in the other dataset
genus_sum_top20$top20 <- 'yes'
genus_sum_top20[genus_sum_top20$rank>20,]$top20 <- 'no'


# plot rank
plot_rank_dia <- ggplot(genus_sum_top20[genus_sum_top20$normalization=='diatoms',], 
                        aes(x = 1, label = rank, y = genus, color = top20)) + 
  geom_text() +
  ggtitle("Rank normalized\nby diatoms") + 
  theme(
    legend.position = "none", 
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "white"),     # Set x-axis labels to white
    axis.title.x = element_text(color = "white")      # Set x-axis title to white
  ) +
  facet_grid(. ~ method) +
  scale_color_manual(values = c("yes" = "black", "no" = "gray"))


############# Plot all panels in a single figure ################


fig6B <- plot_grid(p_classes, plot_rank_dia, p_genus_dia, plot_occupancy, 
                  nrow = 1, 
                  rel_widths=c(0.7,0.5,1,1))

fig6B
#ggsave(fig6B, file='fig6B.pdf')

