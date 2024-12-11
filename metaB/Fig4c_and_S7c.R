
library(mgcv)
library(data.table)
library(tidyverse)
library(ggplot2)
library(vegan)
library(cowplot)
library(scatterpie)

################## Formating the barcode data ################## 

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
rm(v4_processed, v9_processed, process_diatom_data)

######################################################################
############## Physicochemical variables #############################
# original data from https://doi.org/10.1594/PANGAEA.875582
physicochemistry <- read.csv(file='datasets/physicochemistry.tsv', header = T, sep='\t')

##### merge barcode relative abundances with physicochemical variables #####
myabundance$station <- gsub(x=myabundance$station, pattern = '^00', replacement = '')
myabundance$station <- gsub(x=myabundance$station, pattern = '^0', replacement = '')
myabundance <- merge(x=myabundance, y=physicochemistry, by=c('station', 'depth'), all.x = T)
rm(physicochemistry)

##### add the lower pore filter size as an extra variable:
myabundance$lower_size <- sub(myabundance$size_fraction, pattern = '-.*', replacement = '')
myabundance$lower_size <- sub(myabundance$lower_size, pattern = '/.*', replacement = '')
myabundance$lower_size <- as.numeric(as.character(myabundance$lower_size))

############################### NMDS analysis and plot ####################

# define function for NMDS to apply to each barcode method:
calculador <- function(MYMETHOD){  
  mysubset <-  myabundance[myabundance$depth=="SRF",]
  mysubset <-  mysubset[mysubset$method==MYMETHOD,]
  
  mysubset$station <- gsub(mysubset$station, pattern = "^", replacement = "TARA_")
  row.names(mysubset) <- paste(mysubset$station, '_', mysubset$size_fraction, sep='')
  
  vare.mds <- metaMDS(mysubset[,c(8:11)], trymax=10000) # colums 8 to 11 are Araphid_pennate, Mediophyceae, Coscinodiscophyceae, and Raphid_pennate
  print(vare.mds)
  
  #extract NMDS scores (x and y coordinates)
  data.scores <- as.data.frame(scores(vare.mds)$sites) 
  
  #Add columns from your original data (pc) to your new NMDS coordinates data frame. 
  #This will come in handy when you plot your data and want to differentiate groups or treatments:
  data.scores$station <- mysubset$station
  
  data.scores$'Araphid_pennate' <- mysubset$'Araphid_pennate'
  data.scores$'Mediophyceae' <- mysubset$'Mediophyceae'
  data.scores$'Coscinodiscophyceae' <- mysubset$'Coscinodiscophyceae'
  data.scores$'Raphid_pennate' <- mysubset$'Raphid_pennate'
  #data.scores$'unclassified' <- mysubset$'unclassified'
  
  #Using the scores function from vegan to extract the species scores and convert to a data.frame
  species.scores <- as.data.frame(scores(vare.mds, "species"))
  species.scores$species <- rownames(species.scores)
  
  
  data.scores_with_env <- merge(x=data.scores, y=mysubset, all.x=T, by="row.names")
  enviro <- data.scores_with_env %>% select(Row.names, lower_size, Temperature, Ammonium.5m, NO2.5m, NO3.5m, NH4toDIN.5m, Iron.5m, Si, PO4, ChlorophyllA, abslat)
  
  
  enviro <- unique(enviro) 
  rm(data.scores_with_env)
  
  data.envfit <- envfit(vare.mds, enviro, na.rm=T)
  print(data.envfit)
  
  en_coord_cont = as.data.frame(scores(data.envfit, "vectors")) * ordiArrowMul(data.envfit)
  
  
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
  p <- ggplot() +
    geom_scatterpie(data = data.scores, aes(x = NMDS1, y = NMDS2, group = station), pie_scale = (abs(myMIN) + abs(myMAX))/2,
                    cols=c( 'Araphid_pennate', 'Mediophyceae', 'Coscinodiscophyceae', 'Raphid_pennate'), alpha=.85, color=NA) + #unclassified
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
    
    theme(plot.title = element_text(size = 20, face = "bold")) + theme_bw() + ggtitle(MYMETHOD)
  
  p
  
  library(stringr)
  # Extract size fraction
  data.scores$size_fraction <- str_extract(rownames(data.scores), "(?<=_)[0-9./-]+$")
  
  # Normalize size fractions by replacing '/' with 'or'
  data.scores$size_fraction <- str_replace_all(data.scores$size_fraction, "/", "or")
  
  # Plotting
  p_size <- ggplot() +
    # Points colored by normalized size fraction
    geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, color = size_fraction), size = 4, alpha = 0.85) +
    
    # Add species labels
    geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species), color = "blue", show.legend = FALSE) +
    
    # Plot significant vectors (p < 0.05) in red
    geom_segment(data = significant_en_coord, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                 arrow = arrow(length = unit(0.5, "cm")), colour = "red", inherit.aes = FALSE) +
    geom_text(data = significant_en_coord, aes(x = NMDS1, y = NMDS2, label = row.names(significant_en_coord)), color = "red", size = 5) +
    
    # Plot marginally significant vectors (0.05 < p < 0.1) in orange
    geom_segment(data = marginal_en_coord, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                 arrow = arrow(length = unit(0.5, "cm")), colour = "orange", inherit.aes = FALSE) +
    geom_text(data = marginal_en_coord, aes(x = NMDS1, y = NMDS2, label = row.names(marginal_en_coord)), color = "orange", size = 5) +
    
    # Add color legend and theme
    labs(color = "Size Fraction") +
    theme(plot.title = element_text(size = 20, face = "bold")) +
    theme_bw() +
    ggtitle(MYMETHOD)
  
  # Print the plot
  print(p_size)
  
  library(cowplot)
  myplot <- plot_grid(p, p_size)
  
  return(myplot)
  rm(data.scores, species.scores, en_coord_cont, en_coord_cat, data.envfit, myMIN, myMAX, mysubset)
}

plotV9_Fig4c <- calculador("V9")
plotV4_FigS7c <- calculador("V4")


