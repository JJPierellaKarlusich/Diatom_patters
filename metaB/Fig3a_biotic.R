
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
rm("v4", "v4_processed", "p4", "p9", "process_rel_abundance",
    "v4_rel_abund", "v9", "v9_processed", "v9_rel_abund")

# Transform the dataframe to wide format
silicifiers <- silicifiers %>%
  select(-ocean_region) %>% # Remove the 'ocean_region' column
  pivot_wider(names_from = taxon, values_from = percEuk, values_fill = list(percEuk = 0))

silicifiers$station <- gsub(silicifiers$station, pattern = '^00', replacement = '')
silicifiers$station <- gsub(silicifiers$station, pattern = '^0', replacement = '')

#  merge datasets
data_df <- merge(x=silicifiers, y=optical_data, by=c('station', 'depth'), all.x=T)
rm(silicifiers, optical_data)

###############################################################################################
################ PLS analysis #################################################################
###############################################################################################

data_df <- data_df %>% select(station, depth, size_fraction, method, Bacillariophyta, Prochlorococcus, Synechococcus, Copepoda, Centroheliozoa, Choanoflagellida, Chrysophyceae, Dictyochophyceae, Nassellaria, Phaeodarea, Spumellaria, Rhizaria) #Porifera
# mean of repetitions
data_df <- aggregate(data=data_df, FUN=mean, cbind(Bacillariophyta, Prochlorococcus, Synechococcus, Copepoda, Centroheliozoa, Choanoflagellida, Chrysophyceae, Dictyochophyceae, Nassellaria, Phaeodarea, Spumellaria, Rhizaria) ~ station + depth + size_fraction + method, drop_na=F)
data_df <- data_df %>% select(Bacillariophyta, Prochlorococcus, Synechococcus, Copepoda, Centroheliozoa, Choanoflagellida, Chrysophyceae, Dictyochophyceae, Nassellaria, Phaeodarea, Spumellaria, Rhizaria) #Porifera

# Clean the data to remove rows with NAs
clean_data_df <- na.omit(data_df)
#rm(data_df)
# clean_data_df <- decostand(clean_data_df, method = "hellinger")

# Specify the predictors and the response
predictors <- as.matrix(clean_data_df[, -1]) # All columns except the first as predictors
response <- as.matrix(clean_data_df[, 1]) # The first column as the response

# Perform PLS regression
pls_model <- plsreg1(predictors, response, comps = 2)
rm(predictors, response)

# Summary of the model
summary(pls_model)

# plot the circle of correlations
Fig3aa <- as.grob(~plot(pls_model, comps = 1:2))


###############

# Select the first column (Bacillariophyta) and the other columns
all_vars <- data_df[, c("Bacillariophyta", "Prochlorococcus", "Synechococcus", "Copepoda", "Centroheliozoa", "Choanoflagellida", "Chrysophyceae", "Dictyochophyceae", "Nassellaria", "Phaeodarea", "Spumellaria", "Rhizaria")]

# Calculate Spearman correlations and p-values for Bacillariophyta against all other variables using Hmisc::rcorr
cor_results_all <- rcorr(as.matrix(all_vars), type = "spearman")
rm(all_vars)

# Extract correlation coefficients and p-values for Bacillariophyta against the other variables
correlations_all <- cor_results_all$r[1, -1] # First row, excluding the first column (self-correlation)
p_values_all <- cor_results_all$P[1, -1] # First row, excluding the first column (self-correlation)

# Melt the correlation and p-value matrices for ggplot2
melted_correlations_all <- data.frame(Var2 = names(correlations_all), value = correlations_all, p_value = p_values_all)
rm(correlations_all, p_values_all)

# ZooScan: Copepoda, Rhizaria
# FC: Synechococcus, Prochlorococcus

melted_correlations_all$Var2 <- factor(melted_correlations_all$Var2, levels=c("Prochlorococcus", "Synechococcus", "Copepoda", "Centroheliozoa", "Choanoflagellida", "Chrysophyceae", "Dictyochophyceae", "Nassellaria", "Phaeodarea", "Spumellaria", "Rhizaria"))


# Plot the correlations using ggplot2 with a cross for p > 0.05
Fig3ab <- ggplot(melted_correlations_all, aes(x = "Diatom   \nabundance", y = Var2, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name = "Spearman's ρ") +
  theme_minimal() +
  labs(title = "Biotic factors",
       x = "",
       y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = subset(melted_correlations_all, p_value > 0.05), aes(label = "×"), color = "black", size = 5)



#################

Fig3a <- plot_grid(Fig3aa, Fig3ab, ncol=1)
Fig3a
