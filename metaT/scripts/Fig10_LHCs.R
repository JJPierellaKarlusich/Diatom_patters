library(ggplot2)
library(cowplot)
library('data.table')
library(tidyverse)
library(dplyr)
library(vegan)
library(plsdepot)
library(ggplotify)
library(Hmisc)
library(GGally)


# Upload the list of diatom sequences coding for LHCs (PF00504) in version 1.5 of the Marine Atlas of Tara Oceans Unigenes
# The functional classification into the main subfamilies (LHCf, LHCq, LHCr, LHCx and LHCz) was performed by phylogenetic placement 
# of the translated sequences on the reference phylogeny described in Kumazawa et al. 2022 Physiologia Plantarum 174, e13598
lhc <- fread(file='../datasets/LHC/LHC.diatoms.seqs.function.tsv')
lhc$unigene <- sub("MATOU-v1.5.", "", lhc$unigene)
lhc$seq.text <- NULL

# Abundance tables
  # upload metatranscriptomic (metaT) abundance table for diatom unigenes 
  # this table was generated using the scripts in the directory data_generation/
metaT <- fread(file="../datasets/Diatoms/Bacillariophyta.MATOU-v1.5.metaT.occurrences.tsv.gz")
  # select the LHC transcripts
metaT <- metaT[metaT$geneid %in% lhc$unigene,]

  # upload metagenomic (metaG) abundance table for diatom unigenes 
  # this table was generated using the scripts in the directory data_generation/
metaG <- fread(file="../datasets/Diatoms/Bacillariophyta.MATOU-v1.5.metaG.occurrences.tsv.gz")
  # select the LHC genes
metaG <- metaG[metaG$geneid %in% lhc$unigene,]

  # merge metaG and metaT
mydata <- rbind(metaG, metaT)
rm(metaG, metaT)

# select LHC unigenes from the abundance table and merge both dataframes:
mydata$geneid <- as.character(mydata$geneid)
mydata <- merge(x=lhc, y=mydata, by.x='unigene', by.y='geneid', all.y=T)

# formating
mydata$station <- gsub(mydata$sample, pattern="[A-Z].*", replacement="")
mydata$depth <- gsub(mydata$sample, pattern="^[0-9]*", replacement="")
mydata$depth <- gsub(mydata$depth, pattern="[0-9].*", replacement="")
mydata$iteration <- gsub(mydata$sample, pattern="^[0-9]*[A-Z]*", replacement="")
mydata$iteration <- gsub(mydata$iteration, pattern="[A-Z].*", replacement="")
mydata$seqcode <- gsub(mydata$sample, pattern=".*[A-Z]", replacement="")

#selecting the size fraction
mydata$filter <- gsub(mydata$sample, pattern="[0-9]*$", replacement="")
mydata$filter <- gsub(mydata$filter, pattern="^.*[0-9]", replacement="")

mydata$filter <- as.character(mydata$filter)
mydata[mydata$filter=="GGKK",]$filter <- "0.8-3"
mydata[mydata$filter=="GGMM",]$filter <- "0.8-5/2000" #0.8-5
mydata[mydata$filter=="GGQQ",]$filter <- "0.8-5/2000" #"0.8-20"
mydata[mydata$filter=="GGZZ",]$filter <- "0.8-5/2000" #"0.8-2000"
mydata[mydata$filter=="MMQQ",]$filter <- "3/5-20"
mydata[mydata$filter=="QQSS",]$filter <- "20-180"
mydata[mydata$filter=="QQRR",]$filter <- "20-180" #"20-200"
mydata[mydata$filter=="SSUU",]$filter <- "180-2000"
mydata[mydata$filter=="KKQQ",]$filter <- "3/5-20" #"3-20"
mydata[mydata$filter=="KKZZ",]$filter <- "3/5-20" #"3-2000"
mydata <- mydata[mydata$filter!='0.8-3',]

# codes to real info:
mydata[mydata$seqcode=="11",]$seqcode <- "Meta-DNA"
mydata[mydata$seqcode=="12",]$seqcode <- "Meta-DNA_WGA"
mydata[mydata$seqcode=="13",]$seqcode <- "Meta-cDNA_R"
mydata[mydata$seqcode=="14",]$seqcode <- "Meta-cDNA_T"
mydata[mydata$seqcode=="15",]$seqcode <- "Meta-RNA"

mydata$meta <- NA
mydata$meta <- as.character(mydata$meta)
mydata[mydata$seqcode=='Meta-DNA' | mydata$seqcode=='Meta-DNA_WGA',]$meta <- 'metaG'
mydata[mydata$seqcode=='Meta-cDNA_R' | mydata$seqcode=='Meta-cDNA_T'| mydata$seqcode=='Meta-RNA',]$meta <- 'metaT'

#    mydata <- mydata[mydata$seqcode!='Meta-DNA_WGA',]
#    mydata <- mydata[mydata$seqcode!='Meta-RNA',]

# keep only with epipelagic:
mydata <- mydata[mydata$depth=="SUR" | mydata$depth=="DCM",] # | mydata$depth=="MXL"

# sum LHC type per sample:
mydata <- aggregate(data=mydata, FUN=sum, value ~ LHCtype + samplename + station + depth + filter + meta)

# sum of equivalent filters:
mydata <- aggregate(data=mydata, FUN=sum, value ~ LHCtype + station + depth + filter + meta)

# normalize by % LHC
tmp <- aggregate(data=mydata, FUN=sum, value ~ station + depth + filter + meta)
mydata <- merge(x=mydata, y=tmp, by=c("station","depth","filter", 'meta'), all=T)
rm(tmp)
mydata$totalsLHC <- mydata$value.y
mydata$value.y <- NULL
mydata$perc_among_LHC <-100 * mydata$value.x / mydata$totalsLHC

######### Plot size fractions ######################

mydata$filter <- factor(mydata$filter, levels=c('180-2000', '20-180','3/5-20',  '0.8-5/2000'))

size_plot <- ggplot(data=mydata[mydata$LHCtype!='LHCr9Homolog' & mydata$meta=='metaT' & mydata$station!=168,], 
            aes(x = filter, y = perc_among_LHC)) +
  geom_boxplot(outlier.size = 0.5, outlier.shape = 3) + 
  facet_grid(meta ~ LHCtype, scales = 'free') + 
  ylab('% diatom LHC metatranscriptomic reads') + 
  xlab('size fraction (µm)') + 
  coord_flip() 
size_plot


########## Plot SRF vs DCM #########################

dcm <- mydata[mydata$depth=='DCM',]
surf <- mydata[mydata$depth=='SUR',]
surfvsdcm <- merge(x=surf, y=dcm, by=c("station", "filter", "meta", "LHCtype"), all=F)
rm(surf, dcm)
surfvsdcm$filter <- factor(surfvsdcm$filter, levels=c('180-2000', '20-180', '3/5-20', '0.8-5/2000'))
surfvsdcm$log2_SRFtoDCM <- log2(surfvsdcm$perc_among_LHC.x/surfvsdcm$perc_among_LHC.y)
#p <- ggplot(data=surfvsdcm[surfvsdcm$LHCtype!='LHCr9Homolog' & surfvsdcm$meta=='metaT',], aes(x = filter, y = log2_SRFtoDCM, color=filter)) +
#  geom_boxplot(outlier.size = 0.5, outlier.shape = 3) + facet_grid(LHCtype ~ meta) + geom_abline(slope = 0, intercept = 0) + coord_flip()
#p

depth_plot <- ggplot(data=surfvsdcm[surfvsdcm$LHCtype!='LHCr9Homolog' & surfvsdcm$meta=='metaT',], aes(y = LHCtype, x = log2_SRFtoDCM)) +
  geom_boxplot(outlier.size = 0.5, outlier.shape = 3) + 
  facet_grid(meta ~ .) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlab('log2(Surface/DCM)') + ylab('LHC subfamily')
depth_plot

############################################################################

######################################################################
############## PLS LHC types VS physicochemistry ##
######################################################################


# format data:

# original data from https://doi.org/10.1594/PANGAEA.875582
physicochemistry <- read.csv(file='contextual_data/physicochemistry.for.metaT.tsv', header = T, sep='\t')
latlong <- read.csv(file='contextual_data/station_Lat_Long_uniq.withTaraPrefix.tsv', header = T, sep='\t')
physicochemistry <- merge(x=physicochemistry, y=latlong, by='station', all.x = T)
rm(latlong)
physicochemistry$abslat <- abs(physicochemistry$lat)

mydata <- merge(x=mydata, y=physicochemistry, by=c('station', 'depth'), all.x = T)
rm(physicochemistry)

mydata <- mydata[mydata$LHCtype!='LHCr9Homolog',]

mydata$LHCtype <- paste(mydata$meta, mydata$LHCtype, sep='_')
mydata$meta <- NULL

# Convert to wide format
mydata_wide <- mydata %>%
  dplyr::select(-value.x, -totalsLHC) %>%  # Remove value.x and totalsLHC columns
  pivot_wider(
    names_from = LHCtype,           # Use values in LHCtype column for new column names
    values_from = perc_among_LHC    # Use values in perc_among_LHC as values for the new columns
  )

# Replace NA values with 0 for selected columns
mydata_wide <- mydata_wide %>%
  mutate(across(c(metaG_LHCr, metaG_LHCq, metaG_LHCz, metaG_LHCf, metaG_LHCx, metaT_LHCr, metaT_LHCq, metaT_LHCz, metaT_LHCf, metaT_LHCx), ~ replace_na(.x, 0)))


# Create a condition to identify the rows where all the specified columns are equal to 0
condition <- mydata_wide$metaG_LHCr == 0 &
  mydata_wide$metaG_LHCq == 0 &
  mydata_wide$metaG_LHCz == 0 &
  mydata_wide$metaG_LHCf == 0 &
  mydata_wide$metaG_LHCx == 0

# Apply the condition to mutate the relevant columns to NA
mydata_wide <- mydata_wide %>%
  mutate(
    metaG_LHCr = ifelse(condition, NA, metaG_LHCr),
    metaG_LHCq = ifelse(condition, NA, metaG_LHCq),
    metaG_LHCz = ifelse(condition, NA, metaG_LHCz),
    metaG_LHCf = ifelse(condition, NA, metaG_LHCf),
    metaG_LHCx = ifelse(condition, NA, metaG_LHCx)
  )



mydata_wide$size <- gsub(mydata_wide$filter, pattern = '-.*', replacement = '')
mydata_wide$size <- gsub(mydata_wide$size, pattern = '/.*', replacement = '')
mydata_wide$size <- as.numeric(mydata_wide$size)
mydata_wide$layer <- as.numeric(2)
mydata_wide[mydata_wide$depth=='SUR',]$layer <- as.numeric(1)

#LHC_metaT <- mydata_wide[mydata_wide$meta=='metaT',]

clean_data_df <- mydata_wide %>%  dplyr::select(metaG_LHCr, metaG_LHCq, metaG_LHCz, metaG_LHCf, metaG_LHCx, metaT_LHCr, metaT_LHCq, metaT_LHCz, metaT_LHCf, metaT_LHCx, size, layer, Depth.nominal, Temperature, Ammonium.5m, NO2.5m, NO3.5m, Iron.5m, Si, PO4, ChlorophyllA, abslat)

# Clean the data to remove rows with NAs
clean_data_df_noNAs <- na.omit(clean_data_df)

# Define predictors and responses
# 'perc_diatoms' and 'Shannon' are the response variables
# and the rest are the predictors
predictors <- as.matrix(clean_data_df_noNAs[, -(1:10)]) # Excluding the first two columns as they are responses
responses <- as.matrix(clean_data_df_noNAs[, 1:10]) # Only the first two columns are responses
rm(clean_data_df_noNAs)

# Perform PLS regression
# NOTE: The data is scaled to standardized values (mean=0, variance=1) by the plsreg2 function
pls_model_LHC_metaT <- plsreg2(predictors, responses, comps = 2)
rm(predictors, responses)

# Summary of the model
summary(pls_model_LHC_metaT)

# plot the PLS plot
plot_pls_LHC_metaT <- as.grob(~plot(pls_model_LHC_metaT, comps = 1:2))
rm(pls_model_LHC_metaT)


#############################

Fig10efg <- plot_grid(plot_pls_LHC_metaT, size_plot, depth_plot, ncol=1, rel_heights = c(2,1,1))
Fig10efg