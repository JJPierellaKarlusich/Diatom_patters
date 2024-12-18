library(ggplot2)
library(tidyr)
library(tidyverse)
library(data.table)
library(MASS)
library(scales)

# Select here the diatom dataset (Bacillariophyta)
#source('top100_diatoms/metaT_code_with_pfam_description.R')
#rm(depth_abundance, df, OR_abundance, p, p1, p2, p3, p4, size_abundance, df_ranked, pfam_abundance)

# Selecting columns by name
pfam_data_merged <- unique(mydata[c("pfamAcc", "Description")])
rm(mydata)

################## ribosomal proteins:
pfam_data = read.csv(file="contextual_data/PfamA.list",  header =T, sep="\t")
colnames(pfam_data)=c("pfamAcc","Id","Description")
# algunos pfams tienen exactamente la misma descripcion, asi que a esos le agrego su abreviacion para distinguirlos
pfam_data[duplicated(pfam_data$Description),]$Description <- paste(pfam_data[duplicated(pfam_data$Description),]$Id, pfam_data[duplicated(pfam_data$Description),]$Description, sep=" | ")
ribosomal <- pfam_data[grep(pfam_data$Description, pattern="ribosom", ignore.case=T),]
#rm(pfam_data)
ribosomal <- ribosomal[grep(ribosomal$Description, pattern="Mitochondrial", ignore.case=T, invert = T),]
ribosomal <- ribosomal[grep(ribosomal$Description, pattern="Plastid", ignore.case=T, invert = T),]
ribosomal <- ribosomal[ribosomal$pfamAcc!='PF00861',] #Ribosomal L18 of archaea, bacteria, mitoch. and chloroplast


################# vector with top100 Pfams in diatoms
top100 <- as.character(unique(top100Pfams$pfamAcc))
# Split each element by the slash "/"  (e.g. "PF04110/PF16420/PF02991/PF16558")
top100 <- strsplit(top100, "/")
# flatten the list into a single vector
top100 <- unlist(top100)
# Remove "Ribosomal proteins" from the vector
top100 <- top100[top100 != "Ribosomal proteins"]
# add the list of ribosomal Pfams
top100 <- c(top100, ribosomal$pfamAcc)
top100 <- unique(top100)
#rm(ribosomal)


#Replicating this for all other datasets:
final_data = data.frame(pfamAcc = top100Pfams$pfamAcc,  abundance_perc = top100Pfams$abundance_perc, taxon = "Bacillariophyta")

mytaxons <- c("Chlorarachniophyceae", "Chlorophyta", "Cryptophyta", "Dictyochophyceae", "Dinophyceae", "Euglenida", "Haptophyta", "Ochrophyta", "Pelagophyceae")

for (i in mytaxons){
  print(i)
  
  mydata <- read.csv(paste("~/Documents/Tara/Tara_from_website/metaT/datasets/metaT/", i, ".MATOU-v2.Pfam.metaT.tsv" ,sep=""), sep="\t", header=T)

mydata$station <- gsub(mydata$sample, pattern="[A-Z].*", replacement="")
mydata$depth <- gsub(mydata$sample, pattern="^[0-9]*", replacement="")
mydata$depth <- gsub(mydata$depth, pattern="[0-9].*", replacement="")
mydata$iteration <- gsub(mydata$sample, pattern="^[0-9]*[A-Z]*", replacement="")
mydata$iteration <- gsub(mydata$iteration, pattern="[A-Z].*", replacement="")
mydata$seqcode <- gsub(mydata$sample, pattern=".*[A-Z]", replacement="")

#selecting the size fraction
mydata$filter <- gsub(mydata$sample, pattern="[0-9]*$", replacement="")
mydata$filter <- gsub(mydata$filter, pattern="^.*[0-9]", replacement="")

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

# codes to real info:
#mydata[mydata$seqcode=="11",]$seqcode <- "Meta-DNA"
#mydata[mydata$seqcode=="12",]$seqcode <- "Meta-DNA_WGA"
mydata[mydata$seqcode=="13",]$seqcode <- "Meta-cDNA_R"
mydata[mydata$seqcode=="14",]$seqcode <- "Meta-cDNA_T"
mydata[mydata$seqcode=="15",]$seqcode <- "Meta-RNA"

# merge repetitions and equivalent filters (eg 0.8-5 and 0.8-2000)
mydata <- aggregate(data=mydata, FUN=sum, rpkm ~ pfamAcc + taxon + station + depth + seqcode + filter + meta)

# keep only with epipelagic:
mydata <- mydata[mydata$depth=="SUR" | mydata$depth=="DCM" | mydata$depth=="MXL",]

  
  
  # normalization: percentage of abundance of a pfam in a sample
  mydata_sum <- aggregate(data=mydata, rpkm ~ taxon + station + depth + seqcode + filter + meta, FUN=sum) # returns the sum of all protein abundance entries for a given sample
  mydata <- merge(x=mydata, y=mydata_sum, by=c("taxon", "station", "depth", "seqcode", "filter", "meta"), all=T) # on a ajouté l'information de la somme des abondances pour un lieu à toutes les protéines de ce lieu
  rm(mydata_sum)
  mydata$abundance_perc <-  (mydata$rpkm.x / mydata$rpkm.y ) * 100 # adds a new column
  mydata$rpkm.x <- NULL
  mydata$rpkm.y <- NULL
  
  print("Normalisation done")
  
  TOTAL_FOR_THIS_TAXA <- sum(mydata$abundance_perc)
  
  #Select only the top100 pfams of diatoms in these datasets:
  mydata <- mydata[mydata$pfamAcc %in% top100,]
  
  
  #### MERGE SIMILAR PFAMS HERE!!!!!!

  # HMG (high mobility group) box
  mydata[mydata$pfamAcc=="PF00505" | mydata$pfamAcc=="PF09011",]$pfamAcc <- 'PF00505/PF09011'
  
  # Tetratricopeptide repeats (3 pfams entree los 100-200 top, pero en total son 20)
  Tetratricopeptides <- c("PF00515","PF07719","PF07720","PF07721","PF09976","PF13174","PF13176","PF13181","PF13371","PF13374","PF13424","PF13428","PF13429","PF13431","PF13432","PF13512","PF14559","PF14561","PF14852","PF14853","PF16669","PF16918","PF18028","PF18391","PF18710","PF18833") 
  mydata$pfamAcc[mydata$pfamAcc %in% Tetratricopeptides] <- "PF00515/PF07719/PF07720/PF07721/PF09976/PF13174/PF13176/PF13181/PF13371/PF13374/PF13424/PF13428/PF13429/PF13431/PF13432/PF13512/PF14559/PF14561/PF14852/PF14853/PF16669/PF16918/PF18028/PF18391/PF18710/PF18833"
  rm(Tetratricopeptides)
  
  # Pentatricopeptide repeats (2 pfams entre los 100-200 top, pero estan estos dos mas: "PF12854" "PF13041")
  mydata[mydata$pfamAcc=="PF13812" | mydata$pfamAcc=="PF01535" | mydata$pfamAcc=="PF12854" | mydata$pfamAcc=="PF13041",]$pfamAcc <- 'PF13812/PF01535/PF12854/PF13041'
  
  # S-adenosyl-L-homocysteine hydrolase
#  mydata[mydata$pfamAcc=="PF00670" | mydata$pfamAcc=="PF05221",]$pfamAcc <- 'PF00670/PF05221'
  mydata[mydata$pfamAcc %in% c("PF00670", "PF05221"), "pfamAcc"] <- 'PF00670/PF05221'   # Directly updating the rows that match the condition, so I avoid an error in Euglenids as they do not contain these Pfams

  # Pyridine nucleotide-disulphide oxidoreductase
  mydata[mydata$pfamAcc=="PF02852" | mydata$pfamAcc=="PF07992" | mydata$pfamAcc=="PF00070" | mydata$pfamAcc=="PF13738",]$pfamAcc <- 'PF02852/PF07992/PF00070/PF13738'
  
  # Kelch motif (tal vez este tambien: PF07707)
  mydata[mydata$pfamAcc=="PF01344" | mydata$pfamAcc=="PF13418" | mydata$pfamAcc=="PF13964" | mydata$pfamAcc=="PF07646" | mydata$pfamAcc=="PF13854"  | mydata$pfamAcc=="PF07707",]$pfamAcc <- 'PF01344/PF13418/PF13964/PF07646/PF13854/PF07707'
  
  # Ankyrin repeats
  mydata[mydata$pfamAcc=="PF00023" | mydata$pfamAcc=="PF12796" | mydata$pfamAcc=="PF13606" | mydata$pfamAcc=="PF13637" | mydata$pfamAcc=="PF13857",]$pfamAcc <- "PF00023/PF12796/PF13606/PF13637/PF13857"
  
  # Glyceraldehyde 3-phosphate dehydrogenase (NAD binding domain & C-terminal domain)
  mydata[mydata$pfamAcc=="PF00044" | mydata$pfamAcc=="PF02800",]$pfamAcc <- "PF00044/PF02800"
  
  # EF hands (elongation factors)
  efhands <- pfam_data[grep(pfam_data$Description, pattern="EF[- ]hand", ignore.case=T) ,]$pfamAcc
  mydata$pfamAcc[mydata$pfamAcc %in% efhands] <- "PF00036/PF02761/PF08355/PF08356/PF08726/PF08976/PF09068/PF09069/PF12763/PF13202/PF13405/PF13499/PF13833/PF14658/PF14788/PF17901/PF17958/PF17959"
  rm(efhands)
  
  #Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase
  mydata[mydata$pfamAcc=='PF02518' |  mydata$pfamAcc=='PF13589',]$pfamAcc <- 'PF02518/PF13589'
  
  # Elongation factor Tu (hay otro pfam:PF14578)
  mydata[mydata$pfamAcc=='PF00009' |  mydata$pfamAcc=='PF03143' | mydata$pfamAcc=='PF03144' |  mydata$pfamAcc=='PF14578',]$pfamAcc <- 'PF00009/PF03143/PF03144/PF14578'
  
  # Elongation factor G
  mydata[mydata$pfamAcc=='PF00679' | mydata$pfamAcc=='PF03764' | mydata$pfamAcc=='PF14492' | mydata$pfamAcc=='PF07299',]$pfamAcc <- "PF00679/PF03764/PF14492/PF07299"
  
  # Helicase (son 28 pfams!, pero me parece que solo un par son top200)
  # Puede que haya mas helicases en otros fitoplacton, pero no pongo esos Pfams por ahora
  Helicases <- c("PF00270","PF00772","PF10410","PF03796","PF01935","PF16203","PF03457","PF04408","PF00271","PF13307","PF13625","PF05127","PF11408","PF14214","PF02689","PF08482","PF01516","PF05970","PF00519","PF00910","PF05496","PF12513","PF07057","PF09416","PF13361","PF13538","PF00580","PF01443")
  mydata$pfamAcc[mydata$pfamAcc %in% Helicases] <- "PF00270/PF00772/PF10410/PF03796/PF01935/PF16203/PF03457/PF04408/PF00271/PF13307/PF13625/PF05127/PF11408/PF14214/PF02689/PF08482/PF01516/PF05970/PF00519/PF00910/PF05496/PF12513/PF07057/PF09416/PF13361/PF13538/PF00580/PF01443"
  rm(Helicases)
  
  # Methyltransferase
  Methyltransferases <- c("PF00891","PF08241","PF08242","PF12847","PF13383","PF13489","PF13578","PF13649","PF13679","PF13847")
  mydata$pfamAcc[mydata$pfamAcc %in% Methyltransferases] <- "PF00891/PF08241/PF08242/PF12847/PF13383/PF13489/PF13578/PF13649/PF13679/PF13847"
  rm(Methyltransferases)
  
  # Ribosomal proteins (son 3-4 Pfams entre los top200m, pero en total son 140 Pfams, que incluyen uno que es "Nucleolar pre-ribosomal-associated protein 1")
  mydata$pfamAcc[mydata$pfamAcc %in% ribosomal$pfamAcc] <- "Ribosomal proteins"

  # ABC transporter transmembrane region
  # Puede que haya mas ABC transporter transmembrane region, especialmente en otros fitoplacton, pero no pongo esos Pfams por ahora
  mydata[mydata$pfamAcc=="PF00664" | mydata$pfamAcc=="PF00005",]$pfamAcc <- 'PF00664/PF00005'
  
  # ubiquitin (son 41 pfams en total)
  ubiquitins <- c("PF04110","PF16420","PF02991","PF16558","PF08647","PF10302","PF16191","PF16190","PF09358","PF13764","PF00632","PF09814","PF01398","PF04424","PF01088","PF09070","PF04683","PF11976","PF13881","PF16207","PF13019","PF14732","PF14451","PF08587","PF09288","PF10585","PF16455","PF00240","PF14560","PF14377","PF00443","PF13423","PF10399","PF08694","PF03152","PF10408","PF03671","PF02809","PF00179","PF14533","PF02148")
  mydata$pfamAcc[mydata$pfamAcc %in% ubiquitins] <- "PF04110/PF16420/PF02991/PF16558/PF08647/PF10302/PF16191/PF16190/PF09358/PF13764/PF00632/PF09814/PF01398/PF04424/PF01088/PF09070/PF04683/PF11976/PF13881/PF16207/PF13019/PF14732/PF14451/PF08587/PF09288/PF10585/PF16455/PF00240/PF14560/PF14377/PF00443/PF13423/PF10399/PF08694/PF03152/PF10408/PF03671/PF02809/PF00179/PF14533/PF02148"
  rm(ubiquitins)
  
  # Glycosyl hydrolase family (2 pfams en el top200, pero en total son 54 Pfams!!!)
  Glycosylhydrolases <- c("PF00150", "PF00232", "PF00331", "PF00728", "PF00759", "PF00840", "PF00933", "PF01373", "PF01374", "PF01532", "PF01670", "PF01915", "PF02011", "PF02015", "PF02055", "PF02057", "PF02156", "PF02324", "PF02838", "PF03065", "PF03200", "PF03443", "PF03512", "PF03632", "PF03633", "PF03636", "PF03639", "PF03644", "PF03648", "PF03659", "PF03662", "PF03663", "PF03664", "PF03718", "PF07470", "PF07477", "PF07488", "PF07745", "PF07971", "PF08306", "PF08307", "PF13199", "PF13647", "PF14498", "PF14587", "PF14883", "PF14885", "PF15979", "PF16317", "PF16862", "PF16874", "PF16875", "PF16923", "PF17189", "PF17387", "PF17433", "PF17652", "PF17678", "PF18034")
  mydata$pfamAcc[mydata$pfamAcc %in% Glycosylhydrolases] <- "PF00150/PF00232/PF00331/PF00728/PF00759/PF00840/PF00933/PF01373/PF01374/PF01532/PF01670/PF01915/PF02011/PF02015/PF02055/PF02057/PF02156/PF02324/PF02838/PF03065/PF03200/PF03443/PF03512/PF03632/PF03633/PF03636/PF03639/PF03644/PF03648/PF03659/PF03662/PF03663/PF03664/PF03718/PF07470/PF07477/PF07488/PF07745/PF07971/PF08306/PF08307/PF13199/PF13647/PF14498/PF14587/PF14883/PF14885/PF15979/PF16317/PF16862/PF16874/PF16875/PF16923/PF17189/PF17387/PF17433/PF17652/PF17678/PF18034"
  rm(Glycosylhydrolases)
  
  
  
  # merge the equivalent Pfams
  mydata <- aggregate(mydata, FUN=sum, cbind(abundance_perc) ~ .)
  #################
  
  # Summing all samples
  mydata <- aggregate(data=mydata, FUN=sum, abundance_perc ~ pfamAcc + taxon)#pour une protéine je fais la somme de toutes ces abondances relatives dans tous les lieux
  mydata$abundance_perc <- mydata$abundance_perc * 100 / TOTAL_FOR_THIS_TAXA
  rm(TOTAL_FOR_THIS_TAXA)

  
  mydata = mydata[, c(1,3,2)]
  
  print("Sum and select done")
  
  final_data = rbind(final_data, mydata)
  
  rm(mydata, i)
}
rm(mytaxons, ribosomal)

final_data$taxon = factor(final_data$taxon, levels=c("Chlorarachniophyceae", "Chlorophyta", "Euglenida", "Cryptophyta",  "Dinophyceae", "Haptophyta", "Ochrophyta", "Dictyochophyceae", "Pelagophyceae", "Bacillariophyta"))


final_data = merge(final_data, pfam_data_merged, by = "pfamAcc", all.x =TRUE)
#final_data$pfamAcc=NULL

# order by abundance in diatoms
top100Pfams <- merge(x=top100Pfams, y=pfam_data_merged, by="pfamAcc", all.x=T)
final_data$Description <- factor(final_data$Description, levels=top100Pfams[order(top100Pfams$abundance_perc),]$Description) # top100Pfams[order(top100Pfams$abundance_perc),] is the data frame with columns ordered according to abundance croissante


library("RColorBrewer")
mycolors <- brewer.pal(name="Paired", n = 12)

final_data$Group <- "Phytoplankton"


p_phyto = ggplot(final_data)+
  aes(x=Description,y=abundance_perc,fill=taxon)+
  geom_bar(stat = "identity",position="fill")+
# scale_fill_brewer( palette = "Set1", direction=-1) +
  scale_fill_manual( values=mycolors ) + 
  labs(x = "", y = "")+
  scale_y_continuous(labels=scales::percent,expand = c(0,0))+
  coord_flip()+
  theme(axis.text.x = element_text(size=10, vjust=1, hjust=0.5),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(.3,"cm"),
        axis.text.y = element_blank(),
        legend.position = "none") +
  facet_grid(. ~ Group)


p_phyto



