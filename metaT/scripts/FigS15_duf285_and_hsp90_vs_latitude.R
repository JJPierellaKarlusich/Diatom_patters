
library(ggplot2)
library(cowplot)


################### distribution of diatom Pfams
META <- 'metaG'
TAXON <- 'Bacillariophyta'
source("organizador.R")
metaG <- mydata
rm(mydata)

META <- 'metaT'
TAXON <- 'Bacillariophyta'
source("organizador.R")
metaT <- mydata
rm(mydata, META)

Bacillariophyta <- rbind(metaG, metaT)
rm(metaG, metaT)
Bacillariophyta$taxon <- TAXON
rm(TAXON)

#######################
META <- 'metaG'
TAXON <- 'Chlorophyta'
source("organizador.R")
metaG <- mydata
rm(mydata)

META <- 'metaT'
TAXON <- 'Chlorophyta'
source("organizador.R")
metaT <- mydata
rm(mydata, META)

Chlorophyta <- rbind(metaG, metaT)
rm(metaG, metaT)
Chlorophyta$taxon <- TAXON
rm(TAXON)

#######################
META <- 'metaG'
TAXON <- 'Dictyochophyceae'
source("organizador.R")
metaG <- mydata
rm(mydata)

META <- 'metaT'
TAXON <- 'Dictyochophyceae'
source("organizador.R")
metaT <- mydata
rm(mydata, META)

Dictyochophyceae <- rbind(metaG, metaT)
rm(metaG, metaT)
Dictyochophyceae$taxon <- TAXON
rm(TAXON)

#######################
META <- 'metaG'
TAXON <- 'Pelagophyceae'
source("organizador.R")
metaG <- mydata
rm(mydata)

META <- 'metaT'
TAXON <- 'Pelagophyceae'
source("organizador.R")
metaT <- mydata
rm(mydata, META)

Pelagophyceae <- rbind(metaG, metaT)
rm(metaG, metaT)
Pelagophyceae$taxon <- TAXON
rm(TAXON)

#######################
META <- 'metaG'
TAXON <- 'Haptophyta'
source("organizador.R")
metaG <- mydata
rm(mydata)

META <- 'metaT'
TAXON <- 'Haptophyta'
source("organizador.R")
metaT <- mydata
rm(mydata, META)

Haptophyta <- rbind(metaG, metaT)
rm(metaG, metaT)
Haptophyta$taxon <- TAXON
rm(TAXON)


#######################
META <- 'metaG'
TAXON <- 'Dinophyceae'
source("organizador.R")
metaG <- mydata
rm(mydata)

META <- 'metaT'
TAXON <- 'Dinophyceae'
source("organizador.R")
metaT <- mydata
rm(mydata, META)

Dinophyceae <- rbind(metaG, metaT)
rm(metaG, metaT)
Dinophyceae$taxon <- TAXON
rm(TAXON)

###################

mydata <- rbind(Dinophyceae, Haptophyta, Pelagophyceae, Dictyochophyceae, Chlorophyta, Bacillariophyta)
rm(Dinophyceae, Haptophyta, Pelagophyceae, Dictyochophyceae, Chlorophyta, Bacillariophyta)

# To avoid "Domian of unknown function" in each row of the graph:
mydata$Description <- gsub(mydata$Description, pattern="^.*DUF", replacement="DUF")
#mydata$Description <- gsub(mydata$Description, pattern="\\)", replacement="")

# the size fraction where the taxon was most prevalent: 0.8-5/2000, 3/5-20, 20-180, and 180-2000 µm for diatoms and dinoflagellates, 0.8-5/2000 µm for the rest.
tmp1 <- mydata[mydata$filter=='0.8-5/2000',]
tmp2 <- mydata[mydata$taxon=='Bacillariophyta' | mydata$taxon=='Dinophyceae',]
mysubset <- unique(rbind(tmp1, tmp2))
rm(tmp1, tmp2)

############## MAPS


# agregar 0 detection
library(dplyr)
mysamples <- mysubset %>% dplyr::select(station, depth, filter, meta, seqcode, taxon, rpkm.y)
mysamples <- unique(mysamples)

PF03382 <- mysubset[mysubset$pfamAcc=='PF03382',]
PF03382 <- merge(x=mysamples, y=PF03382, all.x=T)
PF03382[is.na(PF03382$rpkm.x),]$rpkm.x <- 0 
PF03382[is.na(PF03382$abundance_perc),]$abundance_perc <- 0 
PF03382[is.na(PF03382$pfamAcc),]$pfamAcc <- 'PF03382'
PF03382[is.na(PF03382$name),]$name <- unique(PF03382[!is.na(PF03382$name),]$name)
PF03382[is.na(PF03382$Description),]$Description <- unique(PF03382[!is.na(PF03382$Description),]$Description)

PF00183 <- mysubset[mysubset$pfamAcc=='PF00183',]
PF00183 <- merge(x=mysamples, y=PF00183, all.x=T)
PF00183[is.na(PF00183$rpkm.x),]$rpkm.x <- 0 
PF00183[is.na(PF00183$abundance_perc),]$abundance_perc <- 0 
PF00183[is.na(PF00183$pfamAcc),]$pfamAcc <- 'PF00183'
PF00183[is.na(PF00183$name),]$name <- unique(PF00183[!is.na(PF00183$name),]$name)
PF00183[is.na(PF00183$Description),]$Description <- unique(PF00183[!is.na(PF00183$Description),]$Description)

tmp <- rbind(PF00183, PF03382)
rm(PF00183, PF03382, mysamples)

latlong <- read.csv(file='analysis/station_Lat_Long_uniq.withTaraPrefix.tsv', sep='\t', header = T)
latlong$Station.label <- NULL

tmp <- merge(x=tmp, y=latlong, by='station', all.x=T)
rm(latlong)
tmp$filter <- factor(tmp$filter, levels=c("0.8-5/2000", "3/5-20", "20-180", "180-2000"))

###############

metaG_PF03382 <- ggplot(data=tmp[tmp$meta=='metaG' & tmp$pfamAcc=='PF03382',], aes(x=lat, y=abundance_perc)) + geom_point(shape=1, alpha=0.5) + geom_smooth() + xlim(-80, 80) + facet_grid(taxon ~ meta, scales = 'free') + xlab('Latitude') + ylab('% reads') + ggtitle('PF03382 DUF285') + theme_bw()
metaT_PF03382 <- ggplot(data=tmp[tmp$meta=='metaT' & tmp$pfamAcc=='PF03382',], aes(x=lat, y=abundance_perc)) + geom_point(shape=1, alpha=0.5) + geom_smooth() + xlim(-80, 80) + facet_grid(taxon ~ meta, scales = 'free') + xlab('Latitude') + ylab('% reads') + ggtitle('PF03382 DUF285') + theme_bw()  
metaG_PF00183 <- ggplot(data=tmp[tmp$meta=='metaG' & tmp$pfamAcc=='PF00183',], aes(x=lat, y=abundance_perc)) + geom_point(shape=1, alpha=0.5) + geom_smooth() + xlim(-80, 80) + facet_grid(taxon ~ meta, scales = 'free') + xlab('Latitude') + ylab('% reads') + ggtitle('PF00183 HSP90') + theme_bw() 
metaT_PF00183 <- ggplot(data=tmp[tmp$meta=='metaT' & tmp$pfamAcc=='PF00183',], aes(x=lat, y=abundance_perc)) + geom_point(shape=1, alpha=0.5) + geom_smooth() + xlim(-80, 80) + facet_grid(taxon ~ meta, scales = 'free') + xlab('Latitude') + ylab('% reads') + ggtitle('PF00183 HSP90') + theme_bw()  

library(cowplot)
figS18 <-plot_grid(metaG_PF03382, metaT_PF03382, metaG_PF00183, metaT_PF00183, ncol=4)
figS18
#ggsave(p, file='latitudinal_gradient_for_duf285_vs_hsp90.pdf')
