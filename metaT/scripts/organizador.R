library(ggplot2)
library(tidyr)
library(tidyverse)
library(data.table)
library(MASS)
library(scales)
library(cowplot)

####################################Initial data processing###################################

# DEFINE THIS:
#META <- 'metaG'
#TAXON <- 'Bacillariophyta'

# Select here the diatom dataset (Bacillariophyta)
mydata <- fread(paste("../datasets/Pfam_sums/", META, "/", TAXON ,".MATOU-v1.5.Pfam.", META, ".tsv.gz", sep=''), sep="\t", header=T)
#rm(META)

# remove PF02100 (Ornithine decarboxylase antizyme) as it's misaasgined to diatoms
mydata <- mydata[mydata$pfamAcc!="PF02100",]

# selecting only SUR (surface) & DCM (deep chlorophyll max) 
mydata <- rbind(mydata[grep(mydata$sample, pattern="DCM"),] , mydata[grep(mydata$sample, pattern="SUR"),] )
mydata$taxon = NULL
#mydata$meta = NULL

# organization of the info:
mydata$station <- gsub(mydata$sample, pattern="[A-Z].*", replacement="")
mydata$depth <- gsub(mydata$sample, pattern="^[0-9]*", replacement="")
mydata$depth <- gsub(mydata$depth, pattern="[0-9].*", replacement="")
mydata$iteration <- gsub(mydata$sample, pattern="^[0-9]*[A-Z]*", replacement="")
mydata$iteration <- gsub(mydata$iteration, pattern="[A-Z].*", replacement="")
mydata$filter <- gsub(mydata$sample, pattern="[0-9]*$", replacement="")
mydata$filter <- gsub(mydata$filter, pattern="^.*[0-9]", replacement="")
mydata$seqcode <- gsub(mydata$sample, pattern=".*[A-Z]", replacement="")

# sample codes to real info:

library(dplyr)
mydata <- mydata %>%
  mutate(seqcode = case_when(
    seqcode == "11" ~ "Meta-DNA",
    seqcode == "12" ~ "Meta-DNA_WGA",
    seqcode == "13" ~ "Meta-cDNA_R",
    seqcode == "14" ~ "Meta-cDNA_T",
    seqcode == "15" ~ "Meta-RNA",
    TRUE ~ seqcode
  ))

#mydata[mydata$seqcode=="11",]$seqcode <- "Meta-DNA"
#mydata[mydata$seqcode=="12",]$seqcode <- "Meta-DNA_WGA"
#mydata[mydata$seqcode=="13",]$seqcode <- "Meta-cDNA_R"
#mydata[mydata$seqcode=="14",]$seqcode <- "Meta-cDNA_T"
#mydata[mydata$seqcode=="15",]$seqcode <- "Meta-RNA"

#mydata[mydata$seqcode=="41",]$seqcode <- "Tags_DNA-18SV4"
#mydata[mydata$seqcode=="42",]$seqcode <- "Tags_DNA-18SV9"
#mydata[mydata$seqcode=="43",]$seqcode <- "Tags_DNA-16SV1V2V3V4"
#mydata[mydata$seqcode=="52",]$seqcode <- "Tags_RNA-18SV9"
#mydata[mydata$seqcode=="61",]$seqcode <- "Tags_dpoP"
#mydata[mydata$seqcode=="62",]$seqcode <- "Tags_mcpP"
#mydata[mydata$seqcode=="72",]$seqcode <- "Tags_DNA-18SV9_WGA"
#mydata[mydata$seqcode=="73",]$seqcode <- "Tags_DNA-16SV1V2V3V4_WGA"
#mydata[mydata$seqcode=="64",]$seqcode <- "Tags_DNA-CeV_600_Big"
#mydata[mydata$seqcode=="65",]$seqcode <- "Tags_DNA-CeV_520_Medium"
#mydata[mydata$seqcode=="66",]$seqcode <- "Tags_DNA-CeV_450_Small"
#mydata[mydata$seqcode=="16",]$seqcode <- "Depleted_Meta-RNA"

#mydata[mydata$filter=="GGKK",]$filter <- "0.8-3"
mydata[mydata$filter=="GGMM",]$filter <- "0.8-5/2000" #0.8-5
mydata[mydata$filter=="GGQQ",]$filter <- "0.8-5/2000" # "0.8-20" (solo en estaciones 47 y 48 SUR, entonces lo paso a 0.8-inf)
mydata[mydata$filter=="GGZZ",]$filter <- "0.8-5/2000" #"0.8-2000"
mydata[mydata$filter=="MMQQ",]$filter <- "3/5-20" #"5-20"
mydata[mydata$filter=="QQSS",]$filter <- "20-180"
mydata[mydata$filter=="SSUU",]$filter <- "180-2000"
mydata[mydata$filter=="KKQQ",]$filter <- "3/5-20" #"3-20"
mydata[mydata$filter=="KKZZ",]$filter <- "3/5-20" #"3-2000"
mydata[mydata$filter=="QQRR",]$filter <- "20-180" #"20-200" (solo en estacion 51 MXL/DCM, que no tiene 20-180)
#mydata[mydata$filter=="CCKK",]$filter <- "0.22-3" 
#mydata[mydata$filter=="AACC",]$filter <- "0-0.2"
#mydata[mydata$filter=="CCII",]$filter <- "0.2-1.6"
#mydata[mydata$filter=="EEGG",]$filter <- "0.45-0.8" 
#mydata[mydata$filter=="CCEE",]$filter <- "0.2-0.45"
#mydata[mydata$filter=="BBCC",]$filter <- "0.1-0.2"
#mydata[mydata$filter=="IIQQ",]$filter <- "1.6-20"
#mydata[mydata$filter=="GGRR",]$filter <- "0.8-200"


######### merge repetitions and equivalent filters (eg 0.8-5 and 0.8-2000)
mydata <- aggregate(data=mydata, FUN=sum, rpkm ~ pfamAcc + station + depth + filter + meta) #+ seqcode


# normalization: percentage of abundance of a pfam in a sample
mydata_sum <- aggregate(data=mydata, rpkm ~ station + depth + filter + meta, FUN=sum) #+ seqcode 
mydata <- merge(x=mydata, y=mydata_sum, by=c("station", "depth",  "filter", "meta"), all=T) #"seqcode",
rm(mydata_sum)
mydata$abundance_perc <-  (mydata$rpkm.x / mydata$rpkm.y ) * 100 # adds a new column

# Pfam description
pfam_data = read.csv(file="contextual_data/PfamA.list",  header =T, sep="\t")
# some pfams have the same exact description, so I add their abbreviation for distinguising them
pfam_data[duplicated(pfam_data$Description),]$Description <- paste(pfam_data[duplicated(pfam_data$Description),]$Id, pfam_data[duplicated(pfam_data$Description),]$Description, sep=" | ")
colnames(pfam_data)=c("pfamAcc","name","Description")
#pfam_data$name=NULL
mydata = merge(mydata, pfam_data, by = "pfamAcc", all.x =TRUE)
rm(pfam_data)


mydata$taxon <- TAXON




