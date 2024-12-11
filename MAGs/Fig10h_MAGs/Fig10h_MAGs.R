library(phylotools)
library(reshape2)
library(tidyverse)
library(cowplot)


# Upload the supplementary table from Delmont et al 2022 CellGenomics
# We need this table for the lineage and the BUSCO completeness 
mags <- read.csv(file='SuppTables_Delmont_et_al_2022_CellGenomics/Table_S03_statistics_nr_SMAGs_METdb.csv', header = T, sep='\t')

# keep only Tara Oceans MAGs
mags <- mags[grep(mags$Genome_Id.final.names, pattern = 'METDB', invert = T),]

# keep only diatom MAGs
mags <- mags[!is.na(mags$Best_taxonomy_PHYLUM) & 
               mags$Best_taxonomy_PHYLUM=='Bacillariophyta',]

# dataframe with diatom class assignation:
mytax <- mags %>% select(Genome_Id.final.names, Best_taxonomy_PHYLUM, Best_taxonomy_CLASS, Best_taxonomy_ORDER, Best_taxonomy_GENRE)
mytax[is.na(mytax$Best_taxonomy_CLASS),]$Best_taxonomy_CLASS <- 'unassigned'
mytax[is.na(mytax$Best_taxonomy_GENRE),]$Best_taxonomy_GENRE <- 'unassigned'
mytax[mytax$Best_taxonomy_GENRE=='New_Bacillariaceae_02',]$Best_taxonomy_GENRE <- 'unassigned'
mytax[mytax$Best_taxonomy_ORDER=='Fragilariales' & !is.na(mytax$Best_taxonomy_ORDER),]$Best_taxonomy_CLASS <- 'Fragilariophyceae'

# dataframe with BUSCO completeness
busco <- mags %>% select(Genome_Id.final.names, BUSCO_completion)

#######################################################################################
########## Gene counts for DUF285 and LHC subfamilies for each diatom MAG #############
#######################################################################################

# DUF285 gene count per MAG
DUF285 <- read.fasta(file ="HHMer_search_on_MAGs/PF03382/PF03382.vs.SMAGs_v1_concat.fna", clean_name = FALSE)
DUF285$Pfam <- 'PF03382_DUF285'
DUF285$mag <- gsub(DUF285$seq.name, pattern = '_[0-9]*\\.[0-9]*\\.[0-9]', replacement ='')
DUF285$mag <- gsub(DUF285$mag, pattern = '_Gene.*', replacement ='')
DUF285$mag <- gsub(DUF285$mag, pattern = '-', replacement ='_')
DUF285$seq.text <- NULL
DUF285 <- unique(DUF285)
DUF285 <- merge(x=DUF285, y=mags, by.x='mag', by.y='Genome_Id.final.names', all.y=T)
DUF285$mag_counts <- 1
DUF285[is.na(DUF285$Pfam),]$mag_counts <- 0
DUF285$Pfam <- 'PF03382_DUF285'

# LHC gene count per MAG
lhc <- read.fasta(file ="HHMer_search_on_MAGs/PF00504/PF00504.vs.SMAGs_v1_concat.fna", clean_name = FALSE)
lhc$Pfam <- 'PF00504_LHC'
lhc$mag <- gsub(lhc$seq.name, pattern = '_[0-9]*\\.[0-9]*\\.[0-9]', replacement ='')
lhc$mag <- gsub(lhc$mag, pattern = '_Gene.*', replacement ='')
lhc$mag <- gsub(lhc$mag, pattern = '-', replacement ='_')
lhc$mag <- gsub(lhc$mag, pattern = '_extended$', replacement ='')
lhc$seq.text <- NULL
lhc <- unique(lhc)
lhc <- merge(x=lhc, y=mags, by.x='mag', by.y='Genome_Id.final.names', all.y=T)
lhc$mag_counts <- 1
lhc$Pfam <- 'PF00504_LHC'

pfam_counts <- rbind(DUF285, lhc)
rm(DUF285) #lhc
rm(mags)

# Gene count for LHC subfamilies per MAG
# Upload the list of diatom sequences coding for LHCs (PF00504) in Tara Oceans MAGs
# The functional classification into the main subfamilies (LHCf, LHCq, LHCr, LHCx and LHCz) was performed by phylogenetic placement 
# of the translated sequences on the reference phylogeny described in Kumazawa et al. 2022 Physiologia Plantarum 174, e13598
LHCtypes <- read.csv(file='LHC_diatomMAGs_annotation/LHC.diatoms.MAGs.function.tsv', header = T, sep='\t')
lhc <- merge(x=lhc, y=LHCtypes, by='seq.name', all.y=T)
rm(LHCtypes)
lhc$Pfam <-lhc$LHCtype
lhc$LHCtype <- NULL

pfam_counts <- rbind(pfam_counts, lhc)
rm(lhc)

# gene family size in each MAG
pfam_counts <- aggregate(data=pfam_counts, FUN=sum, mag_counts ~ mag + Best_taxonomy_PHYLUM + BUSCO_completion + Pfam)

pfam_counts <- pfam_counts[pfam_counts$Pfam!='LHCr9homolog',]


#pfam_counts <- dcast(pfam_counts, mag + Best_taxonomy_PHYLUM + BUSCO_completion ~ Pfam, value.var="mag_counts")
#pfam_counts$LHCr9homolog <- NULL
#pfam_counts[is.na(pfam_counts)] <- 0

#######################################################################################
########## diatom MAG distribution ##############################################
#######################################################################################

# Upload the distribution of the MAGs
mag_distribution <- read.csv(file="SuppTables_Delmont_et_al_2022_CellGenomics/Vertical_coverage_corrected.csv", sep="\t", header=T)
mag_distribution$AAA_Max_Coverage <- NULL
mag_distribution$Size_01 <- NULL
mag_distribution$Size_02 <- NULL
mag_distribution$Size_03 <- NULL
mag_distribution$Size_04 <- NULL
mag_distribution$Arctic <- NULL
mag_distribution$Atlantic <- NULL
mag_distribution$Indian <- NULL
mag_distribution$Mediterranean <- NULL
mag_distribution$Red.Sea <- NULL
mag_distribution$Pacific <- NULL
mag_distribution$Southern <- NULL
mag_distribution <- mag_distribution[-c(1,2,4),]
colnames(mag_distribution) <- mag_distribution[1,]
mag_distribution <- mag_distribution[-1,]

# keep only diatom MAGs
mag_distribution <- mag_distribution[mag_distribution$AAA_SMAG %in% unique(mytax$Genome_Id.final.names),]

mag_distribution <- as.data.frame(t(mag_distribution))
colnames(mag_distribution) <- mag_distribution[1,]
mag_distribution <- mag_distribution[-1,]

mag_distribution$sample <- row.names(mag_distribution)

mag_distribution[,c(1:(ncol(mag_distribution)-1))] <- sapply(mag_distribution[,c(1:(ncol(mag_distribution)-1))],as.numeric)

mag_distribution$sample <- gsub(mag_distribution$sample, pattern="^X", replacement="")

# organization of the info:
mag_distribution$station <- gsub(mag_distribution$sample, pattern="[A-Z].*", replacement="")
mag_distribution$depth <- gsub(mag_distribution$sample, pattern="^[0-9]*", replacement="")
mag_distribution$depth <- gsub(mag_distribution$depth, pattern="[0-9].*", replacement="")
mag_distribution$iteration <- gsub(mag_distribution$sample, pattern="^[0-9]*[A-Z]*", replacement="")
mag_distribution$iteration <- gsub(mag_distribution$iteration, pattern="[A-Z].*", replacement="")
mag_distribution$filter <- gsub(mag_distribution$sample, pattern="[0-9]*$", replacement="")
mag_distribution$filter <- gsub(mag_distribution$filter, pattern="^.*[0-9]", replacement="")
mag_distribution$seqcode <- gsub(mag_distribution$sample, pattern=".*[A-Z]", replacement="")

mag_distribution[mag_distribution$seqcode=="11",]$seqcode <- "Meta-DNA"

#mag_distribution[mag_distribution$filter=="GGKK",]$filter <- "0.8-3"
mag_distribution[mag_distribution$filter=="GGMM",]$filter <- "0.8-5/2000" #"0.8-5"
#mag_distribution[mag_distribution$filter=="GGQQ",]$filter <- "0.8-5/2000" #"0.8-20"
mag_distribution[mag_distribution$filter=="GGZZ",]$filter <- "0.8-5/2000" #"0.8-2000"
mag_distribution[mag_distribution$filter=="MMQQ",]$filter <- "3/5-20" #"5-20"
mag_distribution[mag_distribution$filter=="QQSS",]$filter <- "20-180"
mag_distribution[mag_distribution$filter=="QQRR",]$filter <- "20-180" #"20-200"
mag_distribution[mag_distribution$filter=="SSUU",]$filter <- "180-2000"
mag_distribution[mag_distribution$filter=="KKQQ",]$filter <- "3/5-20" #"3-20"
#mag_distribution[mag_distribution$filter=="KKZZ",]$filter <- "3/5-20" #"3-2000"
mag_distribution[mag_distribution$filter=="CCKK",]$filter <- "0.2-1.6/3" #"0.2-3" 
mag_distribution[mag_distribution$filter=="CCII",]$filter <- "0.2-1.6/3" #"0.2-1.6"

mag_distribution$sample <- NULL
mag_distribution$iteration <- NULL

mag_distribution <- aggregate(data=mag_distribution, FUN=mean, . ~ station + depth + filter + seqcode)

# add the lat and long
lat_long <- read.csv(file="../../metaT/scripts/contextual_data/station_Lat_Long_uniq.withTaraPrefix.tsv", sep="\t", header = T)
mag_distribution <- merge(x=mag_distribution, y=lat_long, by="station", all.x = T)
rm(lat_long)

mag_distribution$depth <- factor(mag_distribution$depth, levels=c("SUR", "DCM"))
mag_distribution$filter <- factor(mag_distribution$filter, levels=c("0.2-1.6/3", "0.8-5/2000", "3/5-20", "20-180", "180-2000" ))

library(tidyr)                                                                                                                                        
mag_distribution <- gather(data = mag_distribution, key = SMAG_id, value = vertical_coverage_corrected, -c(station, depth, filter, seqcode, lat, long))

#######################################################################################
########## plots ##############################################
#######################################################################################

# only keep MAGs with BUSCO completeness >70%
busco <- busco[busco$BUSCO_completion>=70,]
pfam_counts <- pfam_counts[pfam_counts$mag %in% busco$Genome_Id.final.names,]
mytax <- mytax[mytax$Genome_Id.final.names %in% busco$Genome_Id.final.names,]
mag_distribution <- mag_distribution[mag_distribution$SMAG_id %in% busco$Genome_Id.final.names,]

# same MAG order in all dataframes
orderedlist <- c(
  'TARA_IOS_50_MAG_00115', 'TARA_PSW_86_MAG_00261',
  'TARA_MED_95_MAG_00394', 'TARA_PSW_86_MAG_00236',
  'TARA_AON_82_MAG_00159', 'TARA_IOS_50_MAG_00056',
  'TARA_MED_95_MAG_00467', 'TARA_MED_95_MAG_00424',
  'TARA_PSW_86_MAG_00222', 'TARA_PSW_86_MAG_00256',
  'TARA_MED_95_MAG_00399', 'TARA_PSE_93_MAG_00253',
  'TARA_AOS_82_MAG_00050', 'TARA_AON_82_MAG_00338',
  'TARA_ARC_108_MAG_00209', 'TARA_ARC_108_MAG_00219',
  'TARA_SOC_28_MAG_00031', 'TARA_AOS_82_MAG_00176',
  'TARA_ARC_108_MAG_00189', 'TARA_ARC_108_MAG_00187',
  'TARA_ARC_108_MAG_00212', 'TARA_ARC_108_MAG_00217',
  'TARA_PSE_93_MAG_00171', 'TARA_SOC_28_MAG_00037',
  'TARA_ARC_108_MAG_00122', 'TARA_ARC_108_MAG_00108',
  'TARA_ARC_108_MAG_00138', 'TARA_ARC_108_MAG_00253',
  'TARA_ARC_108_MAG_00137', 'TARA_ARC_108_MAG_00116',
  'TARA_ARC_108_MAG_00267', 'TARA_ARC_108_MAG_00165',
  'TARA_SOC_28_MAG_00049', 'TARA_SOC_28_MAG_00077'
)
mytax$Genome_Id.final.names <- factor(mytax$Genome_Id.final.names, levels=orderedlist)
mag_distribution$SMAG_id <- factor(mag_distribution$SMAG_id, levels=orderedlist)
busco$Genome_Id.final.names <- factor(busco$Genome_Id.final.names, levels=orderedlist)
pfam_counts$mag <- factor(pfam_counts$mag, levels=orderedlist)
rm(orderedlist)

mytax$header <- 'Genus'
p0 <- ggplot(data=mytax, aes(y=Genome_Id.final.names, x=1, label=Best_taxonomy_GENRE, color=Best_taxonomy_CLASS)) + 
  geom_text(size=3) +
  theme(legend.position = 'none',
        axis.title.x = element_text(color = "white"),
        axis.text.x = element_text(color = "white"),
        axis.title.y = element_blank()) +
  facet_grid(. ~ header)
p0

mag_distribution$vertical_coverage_corrected <- as.numeric(as.character(mag_distribution$vertical_coverage_corrected))
mag_distribution <- mag_distribution[mag_distribution$vertical_coverage_corrected!=0,]
mag_distribution$header <- 'Latitude'
p1 <- ggplot(data=mag_distribution, aes(y=SMAG_id, x=lat, size=vertical_coverage_corrected)) + 
  geom_point(shape=1) + 
  scale_size_area() + facet_grid(. ~ header) + theme(legend.position = 'none')
p1 <- p1 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
p1

busco$header <- 'BUSCO'
p2 <- ggplot(data=busco, aes(y=Genome_Id.final.names, x=BUSCO_completion)) + 
  geom_bar(stat="identity", position="dodge") +
  facet_grid(. ~ header) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank()) +
  xlab('% completion')
p2

p3 <- ggplot(data=pfam_counts, aes(y=mag, x=mag_counts)) + 
  geom_point(shape=16) + geom_segment(aes(x=0, xend=mag_counts, y=mag, yend=mag)) +
  facet_grid(. ~ Pfam, scale='free')
p3 <- p3 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  xlab('Gene copies')
p3 

Fig10h <- plot_grid(p0, p1, p2, p3, ncol=4, rel_widths = c(1, 0.9, 0.3,2))
Fig10h

