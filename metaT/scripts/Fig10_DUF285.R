
library(ggplot2)
library(cowplot)
library(vegan)
library(plsdepot)
library(ggplotify)
library(Hmisc)
library(GGally)

### DUF285 data
META <- 'metaT'
TAXON <- 'Bacillariophyta'
source("organizador.R")
metaT <- mydata
rm(mydata)

META <- 'metaG'
source("organizador.R")
metaG <- mydata
rm(mydata)
rm(META, TAXON)

mydata <- rbind(metaG, metaT)
rm(metaG, metaT)

# Add 0 detection for DUF285
library(dplyr)
mysamples <- mydata %>% dplyr::select(station, depth, filter, meta, taxon, rpkm.y)
mysamples <- unique(mysamples)

# DUF285 (PF03382)
DUF285 <- mydata[mydata$pfamAcc=='PF03382',]
rm(mydata)
DUF285 <- merge(x=mysamples, y=DUF285, all.x=T)
rm(mysamples)
DUF285[is.na(DUF285$rpkm.x),]$rpkm.x <- 0 
DUF285[is.na(DUF285$abundance_perc),]$abundance_perc <- 0 
DUF285[is.na(DUF285$pfamAcc),]$pfamAcc <- 'PF03382'
DUF285[is.na(DUF285$name),]$name <- unique(DUF285[!is.na(DUF285$name),]$name)
DUF285[is.na(DUF285$Description),]$Description <- unique(DUF285[!is.na(DUF285$Description),]$Description)

latlong <- read.csv(file='contextual_data/station_Lat_Long_uniq.withTaraPrefix.tsv', sep='\t', header = T)

DUF285 <- merge(x=DUF285, y=latlong, by='station', all.x=T)
rm(latlong)

## plot map for metaT

# World map setup
mapWorld <- borders("world", colour = NA, fill = "gray80")

# Create map plot
map_DUF285 <- ggplot() + mapWorld + theme_bw() +
  geom_point(data = DUF285[DUF285$depth=='SUR' & DUF285$meta=='metaT', ], 
             aes(y = lat, x = long, size = abundance_perc), 
             shape = 16, color = 'blue', alpha = 0.5) +
  scale_radius(name = '% diatom\nmetaT\nreads') +
  facet_grid(meta ~ .) +
  ggtitle('DUF285') +
  xlab('Longitude') + ylab('Latitude')

map_DUF285

######################################################################
############## PLS DUF285 types VS physicochemistry ##
######################################################################


# format data:

# original data from https://doi.org/10.1594/PANGAEA.875582
physicochemistry <- read.csv(file='contextual_data/physicochemistry.for.metaT.tsv', header = T, sep='\t')

DUF285 <- merge(x=DUF285, y=physicochemistry, by=c('station', 'depth'), all.x = T)
rm(physicochemistry)
DUF285$abslat <- abs(DUF285$lat)


# Convert to wide format
mydata_wide <- DUF285 %>%
  dplyr::select(-rpkm.x, -rpkm.y, lat, long, Station.label) %>%  # Remove these columns
  pivot_wider(
    names_from = meta,           # Use values in meta column for new column names
    values_from = abundance_perc    # Use values in abundance_perc as values for the new columns
  )


mydata_wide$size <- gsub(mydata_wide$filter, pattern = '-.*', replacement = '')
mydata_wide$size <- gsub(mydata_wide$size, pattern = '/.*', replacement = '')
mydata_wide$size <- as.numeric(mydata_wide$size)
mydata_wide$size <- log(mydata_wide$size+1)
mydata_wide$layer <- as.numeric(2)
mydata_wide[mydata_wide$depth=='SUR',]$layer <- as.numeric(1)

clean_data_df <- mydata_wide %>% dplyr::select(metaG, metaT, size, layer, Depth.nominal, Temperature, Ammonium.5m, NO2.5m, NO3.5m, Iron.5m, Si, PO4, ChlorophyllA, abslat)
rm(mydata_wide)

# Clean the data to remove rows with NAs
clean_data_df_noNAs <- na.omit(clean_data_df)

# Define predictors and responses
# 'metaG' and 'metaT' are the response variables
# and the rest are the predictors
predictors <- as.matrix(clean_data_df_noNAs[, -(1:2)]) # Excluding the first two columns as they are responses
responses <- as.matrix(clean_data_df_noNAs[, 1:2]) # Only the first two columns are responses
rm(clean_data_df_noNAs)

# Perform PLS regression
# NOTE: The data is scaled to standardized values (mean=0, variance=1) by the plsreg2 function
pls_model_DUF285 <- plsreg2(predictors, responses, comps = 2)
rm(predictors, responses)

# Summary of the model
summary(pls_model_DUF285)

# plot the PLS plot
plot_pls_DUF285 <- as.grob(~plot(pls_model_DUF285, comps = 1:2))
rm(pls_model_DUF285)

############## Plot both the map and the PLS in the same graph

Fig10ab <- plot_grid(plot_pls_DUF285, map_DUF285, ncol=1)


