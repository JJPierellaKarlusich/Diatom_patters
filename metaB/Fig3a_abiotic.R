library(data.table)
library(tidyverse)
library(ggplot2)
library(vegan)
library(cowplot)
library(plsdepot)
library(ggplotify)
library(tidyr)
library(reshape2)
library(Hmisc)
library(GGally)


################################################################################
###### Calculate diatom relative abundances and Shannon diversity indexes ######
################################################################################

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
  #combined_data <- combined_data[combined_data$total_subset != 0,]
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

rm(mysamples, v4, v9, process_diatom_data)

######################################################################
############## Physicochemistry as predictor variables ###############
######################################################################

# original data from https://doi.org/10.1594/PANGAEA.875582
physicochemistry <- read.csv(file='datasets/physicochemistry.tsv', header = T, sep='\t')

v4v9_env <- merge(x=v4v9, y=physicochemistry, by=c('station', 'depth'), all.x = T)
rm(physicochemistry, v4v9)

clean_data_df <- v4v9_env %>% select(percEuk, Shannon, Temperature, NH4toDIN.5m, Ammonium.5m, NO2.5m, NO3.5m, Iron.5m, Si, PO4, ChlorophyllA, abslat)
rm(v4v9_env)

######################################################################
############## PLS diatom abundance and Shannon VS physicochemistry ##
######################################################################

#clean_data_df$percEuk <- decostand(clean_data_df$percEuk, method = "hellinger")
clean_data_df <- clean_data_df %>% select(percEuk, Shannon, Temperature, NH4toDIN.5m, Ammonium.5m, NO2.5m, NO3.5m, Iron.5m, Si, PO4, ChlorophyllA, abslat)

# Clean the data to remove rows with NAs
clean_data_df_noNAs <- na.omit(clean_data_df)

# Define predictors and responses
# 'perc_diatoms' and 'Shannon' are the response variables
# and the rest are the predictors
predictors <- as.matrix(clean_data_df_noNAs[, -(1:2)]) # Excluding the first two columns as they are responses
responses <- as.matrix(clean_data_df_noNAs[, 1:2]) # Only the first two columns are responses
rm(clean_data_df_noNAs)

# Perform PLS regression
# NOTE: The data is scaled to standardized values (mean=0, variance=1) by the plsreg2 function
pls_model_abiotic <- plsreg2(predictors, responses, comps = 2)
rm(predictors, responses)

# Summary of the model
summary(pls_model_abiotic)

# plot the PLS plot
plot_pls_abiotic <- as.grob(~plot(pls_model_abiotic, comps = 1:2))
rm(pls_model_abiotic)

################################################################################################## 
############## Spearman rho correlations diatom abundance and Shannon VS physicochemistry ########
################################################################################################## 

# Calculate Spearman correlations and p-values using Hmisc::rcorr
cor_results <- rcorr(as.matrix(clean_data_df[, c(1, 2, 3:ncol(clean_data_df))]), type = "spearman")

# Extract correlation coefficients and p-values
correlations <- cor_results$r[1:2, 3:ncol(cor_results$r)] # First two rows against the rest
p_values <- cor_results$P[1:2, 3:ncol(cor_results$P)] # First two rows against the rest
rm(cor_results)

# Melt the correlation and p-value matrices for ggplot2
melted_correlations <- melt(correlations)
melted_p_values <- melt(p_values)
rm(correlations, p_values)

# Merge melted correlations and p-values
melted_correlations$p_value <- melted_p_values$value
rm(melted_p_values)

# Plot the correlations using ggplot2 with a cross for p > 0.05
plot_cor_abiotic <- ggplot(melted_correlations, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Spearman's ρ") +
  theme_minimal() +
  labs(title = "Abiotic factors",
       x = "",
       y = "Abiotic factors") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = subset(melted_correlations, p_value > 0.05), aes(label = "×"), color = "black", size = 5)


rm(clean_data_df, melted_correlations)


################# Plot PLS and spearman correlations for physicochemical factors #######

Fig3aa <- plot_grid(plot_pls_abiotic, plot_cor_abiotic, ncol=2)
rm(plot_pls_abiotic, plot_cor_abiotic)
Fig3aa

########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

######################################################################
############## Biotic factors as predictor variables ###############
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

############# silicifiers from 18S rDNA metabarcoding ##############
source('Fig1_and_S3_SILICIFIERS.R')
silicifiers <- rbind(v4_rel_abund, v9_rel_abund)
rm( "generate_treemaps", "mysamples", "p4" , "p9", "process_rel_abundance",
    "process_silicifier_data", "v4", "v4_processed",
    "v4_rel_abund", "v9", "v9_processed", "v9_rel_abund")

# Transform the dataframe to wide format
silicifiers <- silicifiers %>%
  select(-ocean_region) %>% # Remove the 'ocean_region' column
  pivot_wider(names_from = taxon, values_from = percEuk, values_fill = list(percEuk = 0))

silicifiers$station <- gsub(silicifiers$station, pattern = '^00', replacement = '')
silicifiers$station <- gsub(silicifiers$station, pattern = '^0', replacement = '')

silicifiers[,c(5:ncol(silicifiers))] <- decostand(silicifiers[,c(5:ncol(silicifiers))], method = "hellinger")


#####  merge datasets
data_df <- merge(x=silicifiers, y=optical_data, by=c('station', 'depth'), all.x=T)
rm(silicifiers, optical_data)

######################################################################
############## PLS diatom abundance VS biotic factors ################
######################################################################

data_df <- data_df %>% select(station, depth, size_fraction, method, Bacillariophyta, Prochlorococcus, Synechococcus, Copepoda, Centroheliozoa, Choanoflagellida, Chrysophyceae, Dictyochophyceae, Nassellaria, Phaeodarea, Spumellaria, Rhizaria) #Porifera
# mean of repetitions
data_df <- aggregate(data=data_df, FUN=mean, cbind(Bacillariophyta, Prochlorococcus, Synechococcus, Copepoda, Centroheliozoa, Choanoflagellida, Chrysophyceae, Dictyochophyceae, Nassellaria, Phaeodarea, Spumellaria, Rhizaria) ~ station + depth + size_fraction + method, drop_na=F)
data_df <- data_df %>% select(Bacillariophyta, Prochlorococcus, Synechococcus, Copepoda, Centroheliozoa, Choanoflagellida, Chrysophyceae, Dictyochophyceae, Nassellaria, Phaeodarea, Spumellaria, Rhizaria) #Porifera

# Clean the data to remove rows with NAs
clean_data_df <- na.omit(data_df)

# Specify the predictors and the response
predictors <- as.matrix(clean_data_df[, -1]) # All columns except the first as predictors
response <- as.matrix(clean_data_df[, 1]) # The first column as the response
rm(clean_data_df)

# Perform PLS regression
# NOTE: The data is scaled to standardized values (mean=0, variance=1) by the plsreg2 function
pls_model_biotic <- plsreg1(predictors, response, comps = 2)
rm(predictors, response)

# Summary of the model
summary(pls_model_biotic)

# plot the circle of correlations
pls_plot_biotic <- as.grob(~plot(pls_model_biotic, comps = 1:2))
#title("Biotic factors")

################################################################################################## 
############## Spearman rho correlations diatom abundance VS biotic factors ######################
################################################################################################## 

# Select the first column (Bacillariophyta) and the other columns
all_vars <- data_df[, c("Bacillariophyta", "Prochlorococcus", "Synechococcus", "Copepoda", "Centroheliozoa", "Choanoflagellida", "Chrysophyceae", "Dictyochophyceae", "Nassellaria", "Phaeodarea", "Spumellaria", "Rhizaria")]
rm(data_df)

# Calculate Spearman correlations and p-values for Bacillariophyta against all other variables using Hmisc::rcorr
cor_results_all <- rcorr(as.matrix(all_vars), type = "spearman")
rm(all_vars)

# Extract correlation coefficients and p-values for Bacillariophyta against the other variables
correlations_all <- cor_results_all$r[1, -1] # First row, excluding the first column (self-correlation)
p_values_all <- cor_results_all$P[1, -1] # First row, excluding the first column (self-correlation)
rm(cor_results_all)

# Melt the correlation and p-value matrices for ggplot2
melted_correlations_all <- data.frame(Var2 = names(correlations_all), value = correlations_all, p_value = p_values_all)
rm(correlations_all, p_values_all)

# ZooScan: Copepoda, Rhizaria
# FC: Synechococcus, Prochlorococcus

melted_correlations_all$Var2 <- factor(melted_correlations_all$Var2, levels=c("Prochlorococcus", "Synechococcus", "Copepoda", "Centroheliozoa", "Choanoflagellida", "Chrysophyceae", "Dictyochophyceae", "Nassellaria", "Phaeodarea", "Spumellaria", "Rhizaria"))


# Plot the correlations using ggplot2 with a cross for p > 0.05
coef_biotic <- ggplot(melted_correlations_all, aes(x = "percEuk", y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Spearman's ρ") +
  theme_minimal() +
  labs(title = "Biotic factors",
       x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = subset(melted_correlations_all, p_value > 0.05), aes(label = "×"), color = "black", size = 5)


rm(melted_correlations_all)

################# Plot PLS and spearman correlations for biotic factors #######

Fig3ab <- plot_grid(pls_plot_biotic, coef_biotic, nrow=1)
rm(pls_plot_biotic, coef_biotic)
Fig3ab

################# merge panels for plotting Fig 3a ########

Fig3a <- plot_grid(Fig3aa, Fig3ab, nrow=1)
Fig3a
