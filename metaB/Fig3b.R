

library(data.table)
library(tidyverse)
library(ggplot2)
library(vegan)
library(cowplot)
library(scatterpie)


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

v4v9 <- rbind(v4_processed, v9_processed)
rm(v4_processed, v9_processed)

v4v9$station <- gsub(x=v4v9$station, pattern = '^00', replacement = '')
v4v9$station <- gsub(x=v4v9$station, pattern = '^0', replacement = '')

rm(v4,v9, process_diatom_data)


######################################################################
############## Physicochemistry and biotic factors ###############
######################################################################

############# picocyanobacteria from flow cytometry ##############
# table downloaded and simplified from: https://data.mendeley.com/datasets/p9r9wttjkm/2
fc <- read.csv(file='datasets/Optical_datasets/flow_cytometry_picocyanobacteria.tsv', sep='\t', header = T)
fc$unit <- NULL

############# zooplankton from ZooScan imaging - WP2 net, 200 μm ##############
# table exported from https://ecotaxa.obs-vlfr.fr/prj/377 and https://ecotaxa.obs-vlfr.fr/prj/378
zoo <- read.csv(file='datasets/Optical_datasets/Zooplankton_Abundance_Size.txt', sep='\t', header = T)
zoo <- zoo[zoo$Dataset=='wp2',]
zoo$Station <- gsub(zoo$Station, pattern = '^TARA_00', replacement = '')
zoo$Station <- gsub(zoo$Station, pattern = '^TARA_0', replacement = '')
zoo$Station <- gsub(zoo$Station, pattern = '^TARA_', replacement = '')
zoo <- zoo %>% select(Station, Zooplankton, Copepoda, Rhizaria, Cnidaria, Tunicata)

##### 
optical_data <- merge(x=fc, y=zoo, all=T, by.x='station', by.y='Station')
rm(fc, zoo)

############## Physicochemistry as predictor variables ###############
# original data from https://doi.org/10.1594/PANGAEA.875582
physicochemistry <- read.csv(file='datasets/physicochemistry.tsv', header = T, sep='\t')

##### 
env_variables <- merge(x=physicochemistry, y=optical_data, all=T, by=c('station', 'depth'))
rm(physicochemistry, optical_data)

v4v9_env <- merge(x=v4v9, y=env_variables, by=c('station', 'depth'), all.x = T)
rm(env_variables, v4v9)


#####

# mean of repetitions
data_df <- v4v9_env %>%
  group_by(station, depth, size_fraction, method) %>%
  summarize( percEuk = mean( percEuk ), Temperature = mean( Temperature ), NH4toDIN.5m=mean(NH4toDIN.5m),
             Ammonium.5m=mean(Ammonium.5m),NO2.5m=mean(NO2.5m), NO3.5m=mean(NO3.5m), Iron.5m=mean(Iron.5m), 
             Si=mean(Si), PO4=mean(PO4), ChlorophyllA=mean(ChlorophyllA), abslat=mean(abslat),
             Prochlorococcus=mean(Prochlorococcus),
             Synechococcus=mean(Synechococcus),
             Copepoda=mean(Copepoda),
             Rhizaria=mean(Rhizaria)
             )
data_df <- as.data.frame(data_df)



clean_data_df <- reshape2::dcast(data=data_df, 
                                 station + depth + method + Temperature + NH4toDIN.5m + Ammonium.5m + NO2.5m + NO3.5m + Iron.5m + Si + PO4 + ChlorophyllA + abslat + Prochlorococcus + Synechococcus + Copepoda + Rhizaria
                                 ~ size_fraction, value.var = "percEuk")
clean_data_df$station <- paste(clean_data_df$station, clean_data_df$method, sep='_')


#clean_data_df <- na.omit(clean_data_df)
clean_data_df <- clean_data_df[!is.na(clean_data_df$'0.8-5/2000'),]
clean_data_df <- clean_data_df[!is.na(clean_data_df$'3/5-20'),]
clean_data_df <- clean_data_df[!is.na(clean_data_df$'20-180'),]
clean_data_df <- clean_data_df[!is.na(clean_data_df$'180-2000'),]

mysubset <-  clean_data_df
#mysubset <-  mysubset[mysubset$depth=="SRF",]

mysubset <- mysubset %>% select("station", "0.8-5/2000", "3/5-20", "20-180", "180-2000", "Temperature", "NH4toDIN.5m", "Ammonium.5m", "NO2.5m", "NO3.5m", "Iron.5m", "Si", "PO4", "ChlorophyllA", "abslat", "Prochlorococcus", "Synechococcus", "Copepoda", "Rhizaria")

mysubset <- mysubset[!is.na(mysubset$Temperature),]
mysubset <- mysubset[!is.na(mysubset$NH4toDIN.5m),]
mysubset <- mysubset[!is.na(mysubset$Ammonium.5m),]
mysubset <- mysubset[!is.na(mysubset$NO2.5m),]
mysubset <- mysubset[!is.na(mysubset$NO3.5m),]
mysubset <- mysubset[!is.na(mysubset$Iron.5m),]
mysubset <- mysubset[!is.na(mysubset$Si),]
mysubset <- mysubset[!is.na(mysubset$PO4),]
mysubset <- mysubset[!is.na(mysubset$ChlorophyllA),]
mysubset <- mysubset[!is.na(mysubset$abslat),]
mysubset <- mysubset[!is.na(mysubset$Prochlorococcus),]
mysubset <- mysubset[!is.na(mysubset$Synechococcus),]
mysubset <- mysubset[!is.na(mysubset$Copepoda),]
mysubset <- mysubset[!is.na(mysubset$Rhizaria),]



vare.mds <- metaMDS(mysubset[,c(2:5)], trymax=1000000000)
print(vare.mds )

#extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(vare.mds)$sites) 

#Add columns from your original data (pc) to your new NMDS coordinates data frame. 
#This will come in handy when you plot your data and want to differentiate groups or treatments:
data.scores$station <- mysubset$station

data.scores$'0.8-5/2000' <- mysubset$'0.8-5/2000'
data.scores$'3/5-20' <- mysubset$'3/5-20'
data.scores$'20-180' <- mysubset$'20-180'
data.scores$'180-2000' <- mysubset$'180-2000'

#Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores <- as.data.frame(scores(vare.mds, "species"))
species.scores$species <- rownames(species.scores)


data.scores_with_env <- merge(x=data.scores, y=mysubset, all.x=T) #by="station"
enviro <- data.scores_with_env %>% select(station, Temperature, NH4toDIN.5m, Ammonium.5m, NO2.5m, NO3.5m, Iron.5m, Si, PO4, ChlorophyllA, abslat,
                                          Prochlorococcus, Synechococcus, Copepoda, Rhizaria)

enviro <- unique(enviro) 
rm(data.scores_with_env)

data.envfit <- envfit(vare.mds, enviro, na.rm=T)
print(data.envfit)


# Filter vectors by p-value
significant_vectors <- data.envfit$vectors$pvals < 0.05
marginal_vectors <- data.envfit$vectors$pvals >= 0.05 & data.envfit$vectors$pvals < 0.1
# Extract and scale coordinates
en_coord_cont <- as.data.frame(scores(data.envfit, "vectors")) * ordiArrowMul(data.envfit)
# Create data frames for significant and marginal vectors
significant_en_coord <- en_coord_cont[significant_vectors, ]
marginal_en_coord <- en_coord_cont[marginal_vectors, ]

# extract scores for factors (stations)
en_coord_cat = as.data.frame(scores(data.envfit, "factors")) * ordiArrowMul(data.envfit)

myMIN <- min(c(species.scores$NMDS1, species.scores$NMDS2, data.scores$NMDS1, data.scores$NMDS2, -en_coord_cont$NMDS1, -en_coord_cont$NMDS2))
myMAX <- max(c(species.scores$NMDS1, species.scores$NMDS2, data.scores$NMDS1, data.scores$NMDS2, -en_coord_cont$NMDS1, -en_coord_cont$NMDS2))

# Plotting
Fig3b <- ggplot() +
  geom_scatterpie(data = data.scores, aes(x = NMDS1, y = NMDS2, group = station), pie_scale = (abs(myMIN) + abs(myMAX)),
                  cols = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"), alpha = .85, color = NA) +
  labs(fill = "Size") +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species, color = "blue"), show.legend = FALSE) + 
  
  # Plot significant vectors (p < 0.05) in red
  geom_segment(data = significant_en_coord, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, "cm")), colour = "red", inherit.aes = FALSE) +
  geom_text(data = significant_en_coord, aes(x = NMDS1, y = NMDS2, label = row.names(significant_en_coord)), color = "red", size = 5) +
  
  # Plot marginally significant vectors (0.05 < p < 0.1) in orange
  geom_segment(data = marginal_en_coord, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.5, "cm")), colour = "orange", inherit.aes = FALSE) +
  geom_text(data = marginal_en_coord, aes(x = NMDS1, y = NMDS2, label = row.names(marginal_en_coord)), color = "orange", size = 5) +
  
  theme(plot.title = element_text(size = 20, face = "bold")) + theme_bw()

Fig3b