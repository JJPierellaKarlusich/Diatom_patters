library(ggplot2)
library(tidyr)
library(tidyverse)
library(data.table)
library(MASS)
library(scales)
library(cowplot)

####################################Initial data processing###################################


# Select here the diatom dataset (Bacillariophyta)
#mydata <- read.csv(file.choose(), sep="\t", header=T)
mydata <- read.csv("../datasets/Pfam_sums/metaT/Bacillariophyta.MATOU-v1.5.Pfam.metaT.tsv.gz", sep="\t", header=T)

# remove PF02100 (Ornithine decarboxylase antizyme)
mydata <- mydata[mydata$pfamAcc!="PF02100",]

# remove obsolete Pfams sin May 2020
mydata <- mydata[mydata$pfamAcc!="PF00992" & mydata$pfamAcc!="PF03164" & mydata$pfamAcc!="PF03490" &
                   mydata$pfamAcc!="PF03986"  & mydata$pfamAcc!="PF08802"  & mydata$pfamAcc!="PF10381"  & 
                   mydata$pfamAcc!="PF10914"  & mydata$pfamAcc!="PF14878",]

# contaminations
mydata <- mydata[mydata$pfamAcc!="PF00992",] #Troponin (PF00992)
mydata <- mydata[mydata$pfamAcc!="PF01576",] #Myosin tail (PF01576)
mydata <- mydata[mydata$pfamAcc!="PF16077",] #Spaetzle (PF16077)

# plastid encoded
mydata <- mydata[mydata$pfamAcc!="PF17088",] #Uncharacterised protein family: YCF90 (PF17088)

# mitochondrial encoded
mydata <- mydata[mydata$pfamAcc!="PF00510",] # COX3 | Cytochrome c oxidase subunit III

#  plastid and mitochondrial ribosomes
mydata <- mydata[mydata$pfamAcc!="PF00861",] ##Ribosomal L18 of archaea, bacteria, mitoch. and chloroplast
pfam_data = read.csv(file="contextual_data/PfamA.list",  header =T, sep="\t")
colnames(pfam_data)=c("pfamAcc","Id","Description")
# algunos pfams tienen exactamente la misma descripcion, asi que a esos le agrego su abreviacion para distinguirlos
pfam_data[duplicated(pfam_data$Description),]$Description <- paste(pfam_data[duplicated(pfam_data$Description),]$Id, pfam_data[duplicated(pfam_data$Description),]$Description, sep=" | ")
mito_plastid_ribo <- pfam_data[grep(pfam_data$Description, pattern="ribosom", ignore.case=T),]
rm(pfam_data)
mito_plastid_ribo <- unique(mito_plastid_ribo[grep("Mitochondrial|Plastid", mito_plastid_ribo$Description, ignore.case=TRUE),])
mydata <- mydata[!(mydata$pfamAcc %in% mito_plastid_ribo$pfamAcc),]
rm(mito_plastid_ribo)




# selecting only SUR (surface) & DCM (deep chlorophyll max) 
mydata <- rbind(mydata[grep(mydata$sample, pattern="DCM"),] , mydata[grep(mydata$sample, pattern="SUR"),] )

mydata$taxon = NULL

# organization of the info:
mydata$station <- gsub(mydata$sample, pattern="[A-Z].*", replacement="")
mydata$depth <- gsub(mydata$sample, pattern="^[0-9]*", replacement="")
mydata$depth <- gsub(mydata$depth, pattern="[0-9].*", replacement="")
mydata$iteration <- gsub(mydata$sample, pattern="^[0-9]*[A-Z]*", replacement="")
mydata$iteration <- gsub(mydata$iteration, pattern="[A-Z].*", replacement="")
mydata$filter <- gsub(mydata$sample, pattern="[0-9]*$", replacement="")
mydata$filter <- gsub(mydata$filter, pattern="^.*[0-9]", replacement="")
mydata$seqcode <- gsub(mydata$sample, pattern=".*[A-Z]", replacement="")

# codes to real info:
#mydata[mydata$seqcode=="11",]$seqcode <- "Meta-DNA"
#mydata[mydata$seqcode=="12",]$seqcode <- "Meta-DNA_WGA"
mydata[mydata$seqcode=="13",]$seqcode <- "Meta-cDNA_R"
mydata[mydata$seqcode=="14",]$seqcode <- "Meta-cDNA_T"
mydata[mydata$seqcode=="15",]$seqcode <- "Meta-RNA"
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
mydata <- aggregate(data=mydata, FUN=sum, rpkm ~ pfamAcc + station + depth + seqcode + filter + meta)


# normalization: percentage of abundance of a pfam in a sample
mydata_sum <- aggregate(data=mydata, rpkm ~ station + depth + seqcode + filter + meta, FUN=sum) # returns the sum of all protein abundance entries for a given sample
mydata <- merge(x=mydata, y=mydata_sum, by=c("station", "depth", "seqcode", "filter", "meta"), all=T) # on a ajouté l'information de la somme des abondances pour un lieu à toutes les protéines de ce lieu
rm(mydata_sum)
mydata$abundance_perc <-  (mydata$rpkm.x / mydata$rpkm.y ) * 100 # adds a new column
mydata$rpkm.y <- NULL
  

################## entering the protein description for each pfam:
# Downloaded from: https://www.ebi.ac.uk/interpro/entry/pfam/#table
pfam_data = read.csv(file="contextual_data/PfamA.list",  header =T, sep="\t")

colnames(pfam_data)=c("pfamAcc","Id","Description")
# algunos pfams tienen exactamente la misma descripcion, asi que a esos le agrego su abreviacion para distinguirlos
pfam_data[duplicated(pfam_data$Description),]$Description <- paste(pfam_data[duplicated(pfam_data$Description),]$Id, pfam_data[duplicated(pfam_data$Description),]$Description, sep=" | ")

mydata = merge(mydata, pfam_data, by = "pfamAcc", all.x =TRUE)
rm(pfam_data)


# merge similar Pfams:

mydata[mydata$pfamAcc=='PF00313',]$Description <- 'Cold-shock DNA-binding domain'

# HMG (high mobility group) box
#mydata[grep(mydata$Description, pattern="HMG", ignore.case=T),]$Description <- "HMG (high mobility group) box"
mydata[mydata$pfamAcc=="PF00505" | mydata$pfamAcc=="PF09011",]$Description <- 'HMG-box domain'
mydata[mydata$pfamAcc=="PF00505" | mydata$pfamAcc=="PF09011",]$pfamAcc <- 'PF00505/PF09011'

# Tetratricopeptide repeats (3 pfams entree los 100-200 top, pero en total son 20)
mydata[grep(mydata$Description, pattern="Tetratricopeptide", ignore.case=T),]$Description <- "Tetratricopeptide repeat"
mydata[grep(mydata$Description, pattern="Tetratricopeptide", ignore.case=T),]$pfamAcc <- "PF00515/PF07719/PF07720/PF07721/PF09976/PF13174/PF13176/PF13181/PF13371/PF13374/PF13424/PF13428/PF13429/PF13431/PF13432/PF13512/PF14559/PF14561/PF14852/PF14853/PF16669/PF16918/PF18028/PF18391/PF18710/PF18833"

# Pentatricopeptide repeats (2 pfams entre los 100-200 top, pero estan estos dos mas: "PF12854" "PF13041")
mydata[mydata$pfamAcc=="PF13812" | mydata$pfamAcc=="PF01535" | mydata$pfamAcc=="PF12854" | mydata$pfamAcc=="PF13041",]$Description <- "Pentatricopeptide repeat"
mydata[mydata$pfamAcc=="PF13812" | mydata$pfamAcc=="PF01535" | mydata$pfamAcc=="PF12854" | mydata$pfamAcc=="PF13041",]$pfamAcc <- 'PF13812/PF01535/PF12854/PF13041'

# S-adenosyl-L-homocysteine hydrolase
mydata[mydata$pfamAcc=="PF00670" | mydata$pfamAcc=="PF05221",]$Description <- "S-adenosyl-L-homocysteine hydrolase"
mydata[mydata$pfamAcc=="PF00670" | mydata$pfamAcc=="PF05221",]$pfamAcc <- 'PF00670/PF05221'

# Pyridine nucleotide-disulphide oxidoreductase
mydata[mydata$pfamAcc=="PF02852" | mydata$pfamAcc=="PF07992" | mydata$pfamAcc=="PF00070" | mydata$pfamAcc=="PF13738",]$Description <- 'Pyridine nucleotide-disulphide oxidoreductase'
mydata[mydata$pfamAcc=="PF02852" | mydata$pfamAcc=="PF07992" | mydata$pfamAcc=="PF00070" | mydata$pfamAcc=="PF13738",]$pfamAcc <- 'PF02852/PF07992/PF00070/PF13738'

# Kelch motif (tal vez este tambien: PF07707)
mydata[mydata$pfamAcc=="PF01344" | mydata$pfamAcc=="PF13418" | mydata$pfamAcc=="PF13964" | mydata$pfamAcc=="PF07646" | mydata$pfamAcc=="PF13854"  | mydata$pfamAcc=="PF07707",]$Description <- 'Kelch motif'
mydata[mydata$pfamAcc=="PF01344" | mydata$pfamAcc=="PF13418" | mydata$pfamAcc=="PF13964" | mydata$pfamAcc=="PF07646" | mydata$pfamAcc=="PF13854"  | mydata$pfamAcc=="PF07707",]$pfamAcc <- 'PF01344/PF13418/PF13964/PF07646/PF13854/PF07707'

# Ankyrin repeats
#mydata[grep(mydata$Description, pattern="Ankyrin", ignore.case=T),]$Description <- "Ankyrin repeats"
mydata[mydata$pfamAcc=="PF00023" | mydata$pfamAcc=="PF12796" | mydata$pfamAcc=="PF13606" | mydata$pfamAcc=="PF13637" | mydata$pfamAcc=="PF13857",]$Description <- 'Ankyrin repeats'
mydata[mydata$pfamAcc=="PF00023" | mydata$pfamAcc=="PF12796" | mydata$pfamAcc=="PF13606" | mydata$pfamAcc=="PF13637" | mydata$pfamAcc=="PF13857",]$pfamAcc <- "PF00023/PF12796/PF13606/PF13637/PF13857"

# Glyceraldehyde 3-phosphate dehydrogenase (NAD binding domain & C-terminal domain)
#mydata[grep(mydata$Description, pattern="Glyceraldehyde", ignore.case=T),]$Description <- "Glyceraldehyde 3-phosphate dehydrogenase"
mydata[mydata$pfamAcc=="PF00044" | mydata$pfamAcc=="PF02800",]$Description <- "Glyceraldehyde 3-phosphate dehydrogenase"
mydata[mydata$pfamAcc=="PF00044" | mydata$pfamAcc=="PF02800",]$pfamAcc <- "PF00044/PF02800"

# EF hands (elongation factors)
# hay 11 pfams ademas de los 4 top200: "PF00036" "PF08356" "PF08726" "PF08976" "PF12763" "PF13202" "PF13405" "PF13499" "PF13833" "PF14658" "PF14788"
mydata[grep(mydata$Description, pattern="EF[- ]hand", ignore.case=T) ,]$Description <- "EF hand domain"
mydata[grep(mydata$Description, pattern="EF[- ]hand", ignore.case=T) ,]$pfamAcc <- "PF00036/PF02761/PF08355/PF08356/PF08726/PF08976/PF09068/PF09069/PF12763/PF13202/PF13405/PF13499/PF13833/PF14658/PF14788/PF17901/PF17958/PF17959"
#mydata[mydata$pfamAcc=='PF13202' |  mydata$pfamAcc=='PF00036' | mydata$pfamAcc=='PF13499' | mydata$pfamAcc=='PF13833',]$Description <- 'EF-hand'
#mydata[mydata$pfamAcc=='PF13202' |  mydata$pfamAcc=='PF00036' | mydata$pfamAcc=='PF13499' | mydata$pfamAcc=='PF13833',]$pfamAcc <- 'PF13202/PF00036/PF13499/PF13833'

#Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase
mydata[mydata$pfamAcc=='PF02518' |  mydata$pfamAcc=='PF13589',]$Description <- 'Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase'
mydata[mydata$pfamAcc=='PF02518' |  mydata$pfamAcc=='PF13589',]$pfamAcc <- 'PF02518/PF13589'

# Elongation factor Tu (hay otro pfam:PF14578)
mydata[grep(mydata$Description, pattern="Elongation factor Tu", ignore.case=T),]$Description <- "Elongation factor Tu"
mydata[mydata$pfamAcc=='PF00009' |  mydata$pfamAcc=='PF03143' | mydata$pfamAcc=='PF03144' |  mydata$pfamAcc=='PF14578',]$pfamAcc <- 'PF00009/PF03143/PF03144/PF14578'

# Elongation factor G
mydata[grep(mydata$Description, pattern="Elongation factor G", ignore.case=T),]$Description <- "Elongation factor G"
mydata[mydata$pfamAcc=='PF00679' | mydata$pfamAcc=='PF03764' | mydata$pfamAcc=='PF14492' | mydata$pfamAcc=='PF07299',]$pfamAcc <- "PF00679/PF03764/PF14492/PF07299"

# Helicase (son 28 pfams!, pero me parece que solo un par son top200)
# Puede que haya mas helicases en otros fitoplacton, pero no pongo esos Pfams por ahora
mydata[grep(mydata$Description, pattern="Helicase", ignore.case=T),]$Description <- "Helicase domains"
mydata[grep(mydata$Description, pattern="Helicase", ignore.case=T),]$pfamAcc <- "PF00270/PF00772/PF10410/PF03796/PF01935/PF16203/PF03457/PF04408/PF00271/PF13307/PF13625/PF05127/PF11408/PF14214/PF02689/PF08482/PF01516/PF05970/PF00519/PF00910/PF05496/PF12513/PF07057/PF09416/PF13361/PF13538/PF00580/PF01443"

# Methyltransferase
mydata[grep(mydata$Description, pattern="Methyltransferase domain", ignore.case=T),]$Description <- "Methyltransferase domain"
mydata[grep(mydata$Description, pattern="Methyltransferase domain", ignore.case=T),]$pfamAcc <- "PF00891/PF08241/PF08242/PF12847/PF13383/PF13489/PF13578/PF13649/PF13679/PF13847"

# Ribosomal proteins (son 3-4 Pfams entre los top200m, pero en total son 140 Pfams, que incluyen uno que es "Nucleolar pre-ribosomal-associated protein 1")
mydata[grep(mydata$Description, pattern="ribosom", ignore.case=T),]$Description <- "Ribosomal proteins"
mydata[grep(mydata$Description, pattern="ribosom", ignore.case=T),]$pfamAcc <- "Ribosomal proteins"

# ABC transporter transmembrane region
mydata[mydata$pfamAcc=="PF00664" | mydata$pfamAcc=="PF00005",]$Description <- 'ABC transporter'
mydata[mydata$pfamAcc=="PF00664" | mydata$pfamAcc=="PF00005",]$pfamAcc <- 'PF00664/PF00005'

# ubiquitin (son 41 pfams en total)
mydata[grep(mydata$Description, pattern="ubiquitin", ignore.case=T),]$Description <- "Ubiquitin domains"
mydata[grep(mydata$Description, pattern="ubiquitin", ignore.case=T),]$pfamAcc <- "PF04110/PF16420/PF02991/PF16558/PF08647/PF10302/PF16191/PF16190/PF09358/PF13764/PF00632/PF09814/PF01398/PF04424/PF01088/PF09070/PF04683/PF11976/PF13881/PF16207/PF13019/PF14732/PF14451/PF08587/PF09288/PF10585/PF16455/PF00240/PF14560/PF14377/PF00443/PF13423/PF10399/PF08694/PF03152/PF10408/PF03671/PF02809/PF00179/PF14533/PF02148"

# Glycosyl hydrolase family (2 pfams en el top200, pero en total son 54 Pfams!!!)
mydata[grep(mydata$Description, pattern="Glycosyl hydrolase", ignore.case=T),]$Description <- "Glycosyl hydrolase family"
mydata[grep(mydata$Description, pattern="Glycosyl hydrolase", ignore.case=T),]$pfamAcc <- "PF00150/PF00232/PF00331/PF00728/PF00759/PF00840/PF00933/PF01373/PF01374/PF01532/PF01670/PF01915/PF02011/PF02015/PF02055/PF02057/PF02156/PF02324/PF02838/PF03065/PF03200/PF03443/PF03512/PF03632/PF03633/PF03636/PF03639/PF03644/PF03648/PF03659/PF03662/PF03663/PF03664/PF03718/PF07470/PF07477/PF07488/PF07745/PF07971/PF08306/PF08307/PF13199/PF13647/PF14498/PF14587/PF14883/PF14885/PF15979/PF16317/PF16862/PF16874/PF16875/PF16923/PF17189/PF17387/PF17433/PF17652/PF17678/PF18034"
#mydata[mydata$pfamAcc=='PF00704' | mydata$pfamAcc=='PF00840' ,]$Description <- "Glycosyl hydrolase family"
#mydata[mydata$pfamAcc=='PF00704' | mydata$pfamAcc=='PF00840',]$pfamAcc <- "PF00704/PF00840"
# Sobala, Łukasz F. "Evolution and phylogenetic distribution of endo-α-mannosidase." Glycobiology 33.9 (2023): 687-699.
# The GH99 protein family underwent a considerable expansion and diversification in diatoms (Diatomeae), some of which have a high number of hits per species. An extreme example is Nitzschia sp. RCC80 (Marron et al. 2016), in which 21 predicted proteins (clustered at 95% identity) were found. GH99 protein sequences from diatoms form two well-supported major clades (Supplementary Fig. S1).
# The only diatom whose N-glycome was studied, Phaeodactylum tricornutum (Baïet et al. 2011; Xie et al. 2021), happens to have only one GH99 protein (clustering within Diatomeae GH99 major clade 2) which was omitted from the bioinformatic analyses in these glycomic studies. Such an expansion, and noncanonical motifs 195–199 and 222–226 (see Fig. 2 and the associated Zenodo repository), imply functional diversification of diatom GH99 proteins.

#################

# merge the equivalent Pfams
mydata <- aggregate(mydata, FUN=sum, cbind(rpkm.x, abundance_perc) ~ .)


# Most abundant Pfam
mydata_sum_agg <- aggregate(data=mydata, FUN=sum, abundance_perc ~ pfamAcc)#pour une protéine je fais la somme de toutes ces abondances relatives dans tous les lieux
mydata_sum_agg = mydata_sum_agg[order(mydata_sum_agg$abundance_perc, decreasing = TRUE),]

TOTAL_FOR_THIS_TAXA <- sum(mydata_sum_agg$abundance_perc)

# 100 top Pfams
top100Pfams <- head(mydata_sum_agg, n=100) # prend les 100 premières protéines les plus abondantes
rm(mydata_sum_agg)
top100Pfams$abundance_perc <- 100 * top100Pfams$abundance_perc / TOTAL_FOR_THIS_TAXA
rm(TOTAL_FOR_THIS_TAXA)

# order by abundance
top100Pfams$pfamAcc <- factor(top100Pfams$pfamAcc, levels=top100Pfams[order(top100Pfams$abundance_perc),]$pfamAcc) # top100Pfams[order(top100Pfams$abundance_perc),] is the data frame with columns ordered according to abundance croissante

# selecting the top 100 Pfams in mydata, I filter for the top 100 pfams
mydata <- mydata[mydata$pfamAcc %in% top100Pfams$pfamAcc,] #gives back mydata with only the lines where the Description is in top100
#rm(top100Pfams)


#------------------------------------------------------------------------#
#-------------------------------PLOTS------------------------------------#
#------------------------------------------------------------------------#

#------------------- Plot:  Pfam abundance (all)-----------

pfam_abundance <- aggregate(data=mydata, FUN=sum, abundance_perc ~ pfamAcc + Description) #creates a data frame with the summed abundance percentage for a given pfam
  colnames(pfam_abundance)=c("pfamAcc", "pfam", "number_reads")


df <- pfam_abundance
library(dplyr)
# Add a ranking column
df_ranked <- df %>%
  mutate(rank = rank(-number_reads, ties.method = "min"))

pfam_abundance <- merge(x=pfam_abundance, y=top100Pfams, by="pfamAcc", all.x=T)

# order by abundance
pfam_abundance$pfam <- factor(pfam_abundance$pfam, levels=pfam_abundance[order(pfam_abundance$number_reads),]$pfam)

pfam_abundance$Header <- "Abundance"

p_abundance = ggplot(pfam_abundance) +
  aes(x =pfam , y = abundance_perc) +
  geom_bar(stat="identity",position="stack") +
  labs(x = "", y = "") +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_text(size=10, vjust=0.5, hjust=0.5),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0.3,"cm"),
        axis.text.y = element_text(size=10,vjust=0.5, hjust=1),
        axis.title.y = element_text(size=22, margin = margin(t = 0, r = 25, b = 0, l = 0)),
        axis.title.x = element_text(size=15, margin = margin(t = 30, r = , b = 0, l = 0)),)+
  coord_flip() +
  facet_grid(. ~ Header)

p_abundance



#--------Plot: size fractions

size_abundance <- aggregate(data=mydata, FUN=sum, abundance_perc ~ Description + filter) # sum of all the abundance_perc for a given place? no for a given pfam and size


# order by abundance
size_abundance$Description <- factor(size_abundance$Description, levels=pfam_abundance[order(pfam_abundance$number_reads),]$pfam)

#order size fraction
size_abundance$filter <-  factor(size_abundance$filter, levels=c("0.8-5/2000","3/5-20","20-180","180-2000"))

size_abundance = size_abundance[complete.cases(size_abundance),] # remove the NA sizes, those that do not correspond to the above levels 

size_abundance$Header <- "Size"


p_size = ggplot(size_abundance)+
  aes(x=Description,y=abundance_perc,fill=filter)+
  geom_bar(stat = "identity",position="fill")+
  scale_fill_brewer(palette="Paired")+
  labs(x = "", y = "")+
  scale_y_continuous(labels=scales::percent, expand = c(0,0))+
  coord_flip()+
  theme(axis.text.x = element_text(size=10, vjust=1, hjust=0.5),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.3,"cm"),
        axis.text.y = element_blank())+
  theme(legend.position = "none") +
  facet_grid(. ~ Header)

p_size




#--------Plot: OR


oceans = read.table(file="contextual_data/station_ocean.tsv", header=TRUE)
oceans$station_whole = NULL
OR_abundance = merge(mydata, oceans, by = "station", all.x =TRUE)
rm(oceans)
OR_abundance <- aggregate(data=OR_abundance, FUN=sum, abundance_perc ~ Description + ocean)

# order by abundance
OR_abundance$Description <- factor(OR_abundance$Description, levels=pfam_abundance[order(pfam_abundance$number_reads),]$pfam)

# order Ocean Region
OR_abundance$ocean <- factor(OR_abundance$ocean, levels=c("MS", "RS", "IO", "SAO", "SO", "SPO", "NPO", "NAO", "AO")) # order according to the expedition

OR_abundance$Header <- "Ocean"


pOR = ggplot(OR_abundance)+
  aes(x = Description, y = abundance_perc, fill=ocean)+
  geom_bar(stat = "identity",position="fill")+
  scale_fill_manual(breaks = c("MS","RS","IO","SAO", "SO", "SPO","NPO","NAO","AO"),
                    values=c("darkblue","firebrick4","chocolate1","yellowgreen","black","lightblue3",
                             "yellow", "forestgreen","brown1"))+
  labs(x = "", y = "")+
  scale_y_continuous(labels=scales::percent,expand = c(0,0))+
  coord_flip()+
  theme(axis.text.x = element_text(size=10, vjust=1, hjust=0.5),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.3,"cm"),
        axis.text.y = element_blank())+
  theme(legend.position = "none") +
  facet_grid(. ~ Header)

pOR


#####################Assembling the plots############################


p=plot_grid(p_abundance, p_size, pOR, ncol=3, rel_widths = c(3,1,1), align="h")
p
