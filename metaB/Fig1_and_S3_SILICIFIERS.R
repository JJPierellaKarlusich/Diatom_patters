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
################ silicifier groups #################
#######################################################

# Define a function to add the silicifiers main groups
process_silicifier_data <- function(data, samples) {

  # Filter silicifiers 
  si_groups <- c(";Bacillariophyta", "Chrysophyceae", ";Dictyochophyceae", ";Spumellaria",
                    ";Nassellaria", ";Phaeodarea", ";Choanoflagellida", ";Centroheliozoa", ";Porifera")
  combined_data <- data[Reduce(`|`, lapply(si_groups, function(group) grepl(group, data$taxonomy)))]

  
  # Define taxonomy patterns with corresponding taxon labels
  taxon_mapping <- list(
    "Bacillariophyta" = ";Bacillariophyta",
    "Chrysophyceae" = "Chrysophyceae",
    "Dictyochophyceae" = ";Dictyochophyceae",
    "Spumellaria" = ";Spumellaria", # Nassellaria and Spumellaria are orders of Rhizaria;Radiolaria;Polycystinea
    "Nassellaria" = ";Nassellaria",  # Nassellaria and Spumellaria are orders of Rhizaria;Radiolaria;Polycystinea
    "Phaeodarea" = ";Phaeodarea",   # Phaeodaria are now placed in Rhizaria;Cercozoa
    "Choanoflagellida" = ";Choanoflagellida",  # Opisthokonta;Choanoflagellida;Choanoflagellatea
    "Centroheliozoa" = ";Centroheliozoa",  # Hacrobia;Centroheliozoa
    "Porifera" = "Porifera" # Sea sponges are classified as benthos during their adult life. This is because they anchor themselves to the seafloor, which is the benthic layer of the marine ecosystem. During their larval stage, however, they are plankton drifting on ocean currents.
  )
  
  # Assign main taxon groups
  combined_data$taxon <- as.character(NA)
  for (taxon in names(taxon_mapping)) {
    pattern <- taxon_mapping[[taxon]]
    combined_data$taxon[grep(paste(pattern, collapse = "|"), combined_data$taxonomy)] <- taxon
  }
  rm(taxon_mapping)

  return(combined_data)
}


# Process v4 and v9 data
v4_processed <- process_silicifier_data(v4, mysamples)
v9_processed <- process_silicifier_data(v9, mysamples)
rm(process_silicifier_data)

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
  si <- abundance_data[-1,]
  rm(abundance_data)
  
  # Calculate total eukaryotic read abundance
  euk_data <- data %>% select(samples$sample_id_pangaea)
  euk_data <- t(euk_data)
  ReadAbundance <- rowSums(euk_data, na.rm = TRUE)
  ReadAbundance <- data.frame(ReadAbundance_Euks = ReadAbundance)
  
  # Merge with `si` and add sample metadata
  si <- merge(si, ReadAbundance, by="row.names", all=TRUE)
  si <- merge(si, samples, by.x="Row.names", by.y="sample_id_pangaea", all.x=TRUE)
  
  # Long format and calculate mean of technical replicates
  si_long <- gather(si, key = taxon, value = percEuk, -c(Row.names, station, depth, size_fraction, ocean_region, ReadAbundance_Euks, sample_material, event_date, event_latitude, event_longitude, depth_nominal, marine_biome, biogeo_province,lower_size_fraction, upper_size_fraction, depthplot, biomeplot, polar, abs_lat, uniq, complete, coral_station, sizeplot, station_plot))
  
  si_long$percEuk <- as.numeric(as.character(si_long$percEuk))
  si_long$percEuk <- 100 * si_long$percEuk / si_long$ReadAbundance_Euks
  si_long <- aggregate(data=si_long, FUN=mean, percEuk ~ ocean_region + station + depth + size_fraction + taxon)

  # order taxa
  si_long$taxon <- factor(si_long$taxon,
                          levels=c('Porifera', 'Centroheliozoa', "Choanoflagellida", 'Phaeodarea', 'Nassellaria', 'Spumellaria', 'Chrysophyceae', 'Dictyochophyceae', 'Bacillariophyta' ))
  
  # Organize size fractions and mean equivalent size fractions
  si_long$size_fraction <- recode(si_long$size_fraction, "3-20" = "3/5-20", "5-20" = "3/5-20", "0.8-5" = "0.8-5/2000", "0.8->" = "0.8-5/2000", "0.8-20" = "0.8-5/2000")
  si_long$size_fraction <- factor(si_long$size_fraction, levels = c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))
  si_long <-  aggregate(percEuk ~ taxon + station + depth + size_fraction + ocean_region, data = si_long, FUN = mean)
  
   
  return(si_long)
}



# Process v4 and v9 relative abundance
v4_rel_abund <- process_rel_abundance(v4_processed, v4, mysamples)
v4_rel_abund$method <- 'V4'

v9_rel_abund <- process_rel_abundance(v9_processed, v9, mysamples)
v9_rel_abund$method <- 'V9'


############# boxplots of relative abundances of eukaryotic silicifiers
p4 <- ggplot(v4_rel_abund, aes(x=taxon, y=percEuk)) +
  geom_boxplot(outlier.stroke =0.2, outlier.shape=1, lwd=0.2) +
  facet_grid(method ~ size_fraction) + 
  scale_y_continuous(trans = scales::log1p_trans(), breaks=c(0,5,10,20,40,80)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("% eukaryotic reads") +
  ggtitle('Silicifier relative abundances')
p4

p9 <- ggplot(v9_rel_abund, aes(x=taxon, y=percEuk)) +
  geom_boxplot(outlier.stroke =0.2, outlier.shape=1, lwd=0.2) +
  facet_grid(method ~ size_fraction) + 
  scale_y_continuous(trans = scales::log1p_trans(), breaks=c(0,5,10,20,40,80)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("% eukaryotic reads") +
  ggtitle('Silicifier relative abundances')
p9



############# treemaps of silicifiers

library(treemap)

# Define the treemap generation function
generate_treemaps <- function(si_data, method) {
  
  # Assign high-level taxonomic groups to `taxon2`
  si_data$taxon2 <- NA
  si_data$taxon2 <- dplyr::case_when(
    si_data$taxon %in% "Centroheliozoa" ~ "Haptista",
    si_data$taxon %in% c("Porifera", "Choanoflagellida") ~ "Opisthokonta",
    si_data$taxon %in% c("Phaeodarea", "Nassellaria", "Spumellaria") ~ "Rhizaria",
    si_data$taxon %in% c("Chrysophyceae", "Dictyochophyceae", "Bacillariophyta") ~ "Ochrophyta",
    TRUE ~ si_data$taxon2
  )
  
  # Get unique size fractions
  unique_filters <- unique(si_data$size_fraction)
  
  # Loop over each size fraction and create a treemap
  for (filter_value in unique_filters) {
    
    print(paste(method, filter_value))
    
    # Filter data for the current size fraction
    filtered_data <- si_data[si_data$size_fraction == filter_value, ]
    
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

