library(data.table)
library(tidyverse)
library(ggplot2)

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

# Read sample metadata
# original table from https://zenodo.org/records/7229815
mysamples <- read.csv("datasets/context_general.tsv", sep = "\t")
mysamples <- mysamples %>%
  filter(depth %in% c("SRF", "DCM")) %>%
  filter(!size_fraction %in% c("0.22-1.6", "0.22-3", "0.8-3"))

# Define a function to process the data
process_data <- function(combined_data, samples) {

# Select only relevant samples and compute read abundance
samples <- samples[samples$sample_id_pangaea %in% colnames(combined_data),]
combined_data <- combined_data %>%
  select(amplicon, taxonomy, confidence, total, spread, sequence, samples$sample_id_pangaea)

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

return(combined_data)
}

# Process v4 and v9 data
v4 <- process_data(v4, mysamples)
v9 <- process_data(v9, mysamples)
rm(process_data)


#######################################################
################ phytoplankton groups #################
#######################################################

# Define a function to add phytoplankton group for each ASV
process_phyto_data <- function(data, samples) {
  # Select photosynthetic dinoflagellates
  target_genera <- c("Pileidinium", "Posoniella", "Ceratium", "Neoceratium", "Protoceratium",
                     "Ceratocorys", "Pyrodinium", "Alexandrium", "Amylax", "Gonyaulax", 
                     "Lingulodinium", "Coolia", "Ostreopsis", "Pyrocystis", "Fragilidium", 
                     "Pyrophacus", "Testudodinium", "Brachidinium", "Akashiwo", "Amphidinium",
                     "Apicoporus", "Bispinodinium", "Cochlodinium", "Gymnodinium", 
                     "Gyrodinium", "Lepidodinium", "Spiniferodinium", "Togula", "Torodinium", 
                     "Hemidinium", "Karenia", "Karlodinium", "Takayama", "Polykrikos", 
                     "Jadwigia", "Heterodinium", "Woloszynskia", "Blastodinium", "Adenoides", 
                     "Peridinium", "Scrippsiella", "Peridiniopsis", "Amphidiniella", 
                     "Thecadinium", "Ensiculifera", "Azadinium", "Brandtodinium", 
                     "Bysmatrum", "Glenodiniopsis", "Cachonina", "Heterocapsa", "Durinskia", 
                     "Galeidinium", "Kryptoperidinium", "Herdmania", "Calciodinellum", 
                     "Glenodinium", "Leonella", "Parvodinium", "Pentapharsodinium", 
                     "Theleodinium", "Thoracosphaera", "Gloeodinium", "Phytodinium", 
                     "Prorocentrum", "Biecheleria", "Polarella", "Pelagodinium", 
                     "Protodinium", "Borghiella", "Symbiodinium", "Spatulodinium")
  
  dinos <- data[grep(";Dinophyceae", data$taxonomy)]
  dinosPhoto <- dinos[Reduce(`|`, lapply(target_genera, function(genus) grepl(genus, dinos$taxonomy)))]
  
  # Filter phytoplankton groups
  phyto_groups <- c(";Ochrophyta", ";Archaeplastida", ";Chrompodellids", ";Chlorarachniophyceae",
                    ";Euglyphida", ";Haptophyta", ";Cryptophyta")
  phyto_data <- data[Reduce(`|`, lapply(phyto_groups, function(group) grepl(group, data$taxonomy)))]
  
  # Filter out non-photosynthetic taxa
  non_photosynthetic_taxa <- c("MOCH-3", "MOCH-4", 
                               ";Pteridomonas", ";Ciliophrys_infusionum", 
                               ";Pedospumella", ";Spumella", ";Poteriospumella", 
                               ";Oikomonas", ";Chrysophyceae_Clade-F", 
                               ";Chrysophyceae_Clade-G", ";Chrysophyceae_Clade-H", 
                               ";Chrysophyceae_Clade-I")
  phyto_data <- phyto_data[!Reduce(`|`, lapply(non_photosynthetic_taxa, function(pattern) grepl(pattern, phyto_data$taxonomy)))]
  
  # Combine dinoflagellates and phytoplankton data
  combined_data <- rbind(dinosPhoto, phyto_data)
  
  # Define taxonomy patterns with corresponding taxon labels
  taxon_mapping <- list(
    "Dinophyceae" = ";Dinophyceae",
    "Chrompodellids" = ";Chrompodellids",
    "Chlorarachniophyceae" = ";Chlorarachniophyceae",
    "Euglyphida" = ";Euglyphida",
    "Haptophyta" = ";Haptophyta",
    "Cryptophyta" = ";Cryptophyta",
    "other Ochrophyta" = ";Ochrophyta",
    "Bacillariophyta" = c(";Bacillariophyta", ";Pseudosolenia_calcar\\-avis"),
    "Bolidophyceae" = ";Bolidophyceae",
    "Dictyochophyceae" = ";Dictyochophyceae",
    "Pelagophyceae" = ";Pelagophyceae",
    "Chrysophyceae" = c("Chrysophyceae", "Antarctosaccion"),
    "MOCH-1/2" =c("MOCH-1","MOCH-2"), 
    "MOCH-5" = "MOCH-5",
    "Phaeophyceae" = "Phaeophyceae",
    "Pinguiophyceae" = "Pinguiophyceae",
    "Raphidophyceae" = "Raphidophyceae",
    "Archaeplastida" = ";Archaeplastida",
    "Rhodophyta" = ";Rhodophyta",
    "Streptophyta" = ";Streptophyta",
    "Ulvophyceae" = ";Ulvophyceae",
    "Chlorophyceae" = ";Chlorophyceae",
    "Mamiellophyceae" = ";Mamiellophyceae",
    "Pyramimonadophyceae" = ";Pyramimonadophyceae",
    "Trebouxiophyceae" = ";Trebouxiophyceae",
    "Prasinodermophyta" = ";Prasinodermophyta",
    "Prasino-Clade-VIII" = ";Prasino\\-Clade\\-VIII",
    "Chlorodendrophyceae" = ";Chlorodendrophyceae",
    "Chloropicophyceae" = ";Chloropicophyceae",
    "Nephroselmidophyceae" = ";Nephroselmidophyceae",
    "Pedinophyceae" = ";Pedinophyceae",
    "Prasino-Clade-9" = ";Prasino\\-Clade\\-9",
    "unclassified_Chlorophyta" = c(";unclassified_Chlorophyta", ";Chlorophyta_X"),
    "unclassified_Archaeplastida" = ";unclassified_Archaeplastida"
  )
  
  # Assign main taxon groups
  combined_data$taxon <- as.character(NA)
  for (taxon in names(taxon_mapping)) {
    pattern <- taxon_mapping[[taxon]]
    combined_data$taxon[grep(paste(pattern, collapse = "|"), combined_data$taxonomy)] <- taxon
  }
  
  return(combined_data)
}



# Process v4 and v9 data
v4_processed <- process_phyto_data(v4, mysamples)
v9_processed <- process_phyto_data(v9, mysamples)
rm(process_phyto_data)

#############################
#############################

# Define a function to calculate relative abundances per main taxon
process_rel_abundance <- function(combined_data, data, samples) {
  
# only samples present in the abundance matrix
samples <- samples[samples$sample_id_pangaea %in% colnames(combined_data),]

# Calculate relative abundance matrix by taxon
abundance_data <- combined_data %>%
  select(taxon, samples$sample_id_pangaea) %>%
  group_by(taxon) %>%
  summarise(across(everything(), sum, na.rm = TRUE))

abundance_data <- as.data.frame(t(abundance_data))
colnames(abundance_data) <- abundance_data[1,]
phyto <- abundance_data[-1,]
rm(abundance_data)

# Calculate total eukaryotic read abundance
euk_data <- data %>% select(samples$sample_id_pangaea)
euk_data <- t(euk_data)
ReadAbundance <- rowSums(euk_data, na.rm = TRUE)
ReadAbundance <- data.frame(ReadAbundance_Euks = ReadAbundance)

# Merge with `phyto` and add sample metadata
phyto <- merge(phyto, ReadAbundance, by="row.names", all=TRUE)
phyto <- merge(phyto, samples, by.x="Row.names", by.y="sample_id_pangaea", all.x=TRUE)

# Long format and calculate mean of technical replicates
phyto_long <- gather(phyto, key = taxon, value = percEuk, -c(Row.names, station, depth, size_fraction, ocean_region, ReadAbundance_Euks, sample_material, event_date, event_latitude, event_longitude, depth_nominal, marine_biome, biogeo_province,lower_size_fraction, upper_size_fraction, depthplot, sizeplot, station_plot, biomeplot, polar, abs_lat, uniq, complete, coral_station))

phyto_long$percEuk <- as.numeric(as.character(phyto_long$percEuk))
phyto_long$percEuk <- 100 * phyto_long$percEuk / phyto_long$ReadAbundance_Euks
phyto_long <- aggregate(data=phyto_long, FUN=mean, percEuk ~ ocean_region + station + depth + size_fraction + taxon)

# Group green algal taxonomic levels and aggregate data
phyto_long <- phyto_long %>%
  mutate(taxon = case_when(
    taxon %in% c("Pedinophyceae", "Chlorodendrophyceae", "Trebouxiophyceae", "Ulvophyceae", "Chlorophyceae") ~ "core Chlorophyta",
    taxon %in% c("Pyramimonadophyceae", "Nephroselmidophyceae", "Prasino-Clade-VIII", "Prasino-Clade-9") ~ "Other prasinophytes",
    TRUE ~ taxon
  ))
phyto_long <-  aggregate(percEuk ~ taxon + station + depth + size_fraction + ocean_region, data = phyto_long, FUN = sum)

# order taxa
phyto_long$taxon <- factor(phyto_long$taxon, levels=c(
  "Cryptophyta",
  "Haptophyta", 
  "Rhodophyta",
  "Mamiellophyceae", #clade II
  "Other prasinophytes",  #  "Pyramimonadophyceae", #clade I /   #  "Nephroselmidophyceae", #clade III /   #  "Prasino-Clade-VIII", / #  "Prasino-Clade-9", 
  "Chloropicophyceae",   #clade VII 
  "core Chlorophyta", #"Pedinophyceae", "Chlorodendrophyceae", "Trebouxiophyceae", "Ulvophyceae", "Chlorophyceae", 
  "unclassified_Chlorophyta" ,
  "Streptophyta",
  "Prasinodermophyta", #third phylum within green algae
  "unclassified_Archaeplastida",
  "Euglyphida",
  "Chlorarachniophyceae",
  "Chrompodellids",
  "Dinophyceae",
  "other Ochrophyta",
  "Pinguiophyceae",   
  "Phaeophyceae",
  "Raphidophyceae",
  "Chrysophyceae",
  "MOCH-5",
  "Pelagophyceae", 
  "MOCH-1/2",   
  "Dictyochophyceae",
  "Bolidophyceae",
  "Bacillariophyta"
))

# Organize size fractions and mean equivalent size fractions
phyto_long$size_fraction <- recode(phyto_long$size_fraction, "3-20" = "3/5-20", "5-20" = "3/5-20", "0.8-5" = "0.8-5/2000", "0.8->" = "0.8-5/2000", "0.8-20" = "0.8-5/2000")
phyto_long$size_fraction <- factor(phyto_long$size_fraction, levels = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
phyto_long <-  aggregate(percEuk ~ taxon + station + depth + size_fraction + ocean_region, data = phyto_long, FUN = mean)


return(phyto_long)
}



# Process v4 and v9 relative abundance
v4_rel_abund <- process_rel_abundance(v4_processed, v4, mysamples)
v4_rel_abund$method <- 'V4'

v9_rel_abund <- process_rel_abundance(v9_processed, v9, mysamples)
v9_rel_abund$method <- 'V9'
rm(process_rel_abundance)

############# boxplots of relative abundances of eukaryotic phytoplankton
p4 <- ggplot(v4_rel_abund, aes(x=taxon, y=percEuk)) +
  geom_boxplot(outlier.stroke =0.2, outlier.shape=1, lwd=0.2) +
  facet_grid(method ~ size_fraction) + 
  scale_y_continuous(trans = scales::log1p_trans(), breaks=c(0,5,10,20,40,80)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("% eukaryotic reads") +
  ggtitle('Eukaryotic phytoplankton relative abundances')
p4

p9 <- ggplot(v9_rel_abund, aes(x=taxon, y=percEuk)) +
  geom_boxplot(outlier.stroke =0.2, outlier.shape=1, lwd=0.2) +
  facet_grid(method ~ size_fraction) + 
  scale_y_continuous(trans = scales::log1p_trans(), breaks=c(0,5,10,20,40,80)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("% eukaryotic reads") +
  ggtitle('Eukaryotic phytoplankton relative abundances')
p9

############# treemaps of eukaryotic phytoplankton

library(treemap)

# Define the treemap generation function
generate_treemaps <- function(phyto_data, method) {
  
  # Assign high-level taxonomic groups to `taxon2`
  phyto_data$taxon2 <- NA
  phyto_data$taxon2 <- dplyr::case_when(
    phyto_data$taxon %in% c("Cryptophyta", "Haptophyta") ~ "Hacrobia",
    phyto_data$taxon %in% c("Rhodophyta", "Mamiellophyceae", "Other prasinophytes", 
                            "Chloropicophyceae", "core Chlorophyta", "unclassified_Chlorophyta",
                            "Streptophyta", "Prasinodermophyta", "unclassified_Archaeplastida") ~ "Archaeplastida",
    phyto_data$taxon %in% c("Euglyphida", "Chlorarachniophyceae") ~ "Rhizaria",
    phyto_data$taxon %in% c("Chrompodellids", "Dinophyceae") ~ "Alveolates",
    phyto_data$taxon %in% c("other Ochrophyta", "Pinguiophyceae", "Phaeophyceae", 
                            "Raphidophyceae", "Chrysophyceae", "MOCH-5", "Pelagophyceae", 
                            "MOCH-1/2", "Dictyochophyceae", "Bolidophyceae", "Bacillariophyta") ~ "Ochrophyta",
    TRUE ~ phyto_data$taxon2
  )
  
  # Get unique size fractions
  unique_filters <- unique(phyto_data$size_fraction)
  
  # Loop over each size fraction and create a treemap
  for (filter_value in unique_filters) {
    
    print(paste(method, filter_value))
    
    # Filter data for the current size fraction
    filtered_data <- phyto_data[phyto_data$size_fraction == filter_value, ]
    
    # Aggregate data to calculate mean `percEuk` by relevant groups
    aggregated_data <- aggregate(data = filtered_data, percEuk ~ station + depth + size_fraction + taxon + taxon2, FUN = mean)
    
    # Prepare filename and create PDF output
    pdf_filename <- paste0(method, "_", gsub("/", "or", filter_value), ".pdf")
    pdf(pdf_filename)
    
    # Generate the treemap
    treemap(
      aggregated_data,
      index = c("taxon2", "taxon"),
      vSize = "percEuk",
      vColor = "taxon2",
      type = "categorical",
      drop.unused.levels = FALSE,
      fontsize.labels = 12,
      force.print.labels = TRUE,
      lowerbound.cex.labels = 0.6,
      palette = "RdYlBu",
      title.legend = method,
      algorithm = "pivotSize",
      aspRatio = 0.5,
      title = filter_value
    )
    
    # Close PDF file
    dev.off()
  }
}

# save pdf file for each treemap:
generate_treemaps(v4_rel_abund, "V4")
generate_treemaps(v9_rel_abund, "V9")
rm(generate_treemaps)

############# boxplots of ASV richness for the main eukaryotic phytoplankton

# ASV number per taxon
richness_v4 <- aggregate(data=v4_processed, FUN=length,  amplicon ~ taxon)
richness_v4$method <- 'V4'
richness_v9 <- aggregate(data=v9_processed, FUN=length,  amplicon ~ taxon)
richness_v9$method <- 'V9'
richness <- rbind(richness_v4, richness_v9)
rm(richness_v4, richness_v9)

# Group green algal taxonomic levels and aggregate data
richness <- richness %>%
  mutate(taxon = case_when(
    taxon %in% c("Pedinophyceae", "Chlorodendrophyceae", "Trebouxiophyceae", "Ulvophyceae", "Chlorophyceae") ~ "core Chlorophyta",
    taxon %in% c("Pyramimonadophyceae", "Nephroselmidophyceae", "Prasino-Clade-VIII", "Prasino-Clade-9") ~ "Other prasinophytes",
    TRUE ~ taxon
  ))
richness <- aggregate(data=richness, FUN=sum, amplicon ~ taxon + method)

# order taxa
richness$taxon <- factor(richness$taxon, levels=c(
  "Cryptophyta",
  "Haptophyta", 
  "Rhodophyta",
  "Mamiellophyceae", #clade II
  "Other prasinophytes",  #  "Pyramimonadophyceae", #clade I /   #  "Nephroselmidophyceae", #clade III /   #  "Prasino-Clade-VIII", / #  "Prasino-Clade-9", 
  "Chloropicophyceae",   #clade VII 
  "core Chlorophyta", #"Pedinophyceae", "Chlorodendrophyceae", "Trebouxiophyceae", "Ulvophyceae", "Chlorophyceae", 
  "unclassified_Chlorophyta" ,
  "Streptophyta",
  "Prasinodermophyta", #third phylum within green algae
  "unclassified_Archaeplastida",
  "Euglyphida",
  "Chlorarachniophyceae",
  "Chrompodellids",
  "Dinophyceae",
  "other Ochrophyta",
  "Pinguiophyceae",   
  "Phaeophyceae",
  "Raphidophyceae",
  "Chrysophyceae",
  "MOCH-5",
  "Pelagophyceae", 
  "MOCH-1/2",   
  "Dictyochophyceae",
  "Bolidophyceae",
  "Bacillariophyta"
))

p0 <- ggplot(richness, aes(x=factor(taxon), y=amplicon)) + 
  geom_point() + 
  geom_segment(aes(x=factor(taxon), xend=factor(taxon), y=0, yend=amplicon)) + 
  facet_grid(. ~ method) +
  coord_flip() + 
  theme_bw() + 
  ylab("# ASVs") + xlab("taxon") +
  ggtitle('Eukaryotic phytoplankton richness')

p0






