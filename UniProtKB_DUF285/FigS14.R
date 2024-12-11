
library('data.table')

# Table from: https://www.uniprot.org/uniprotkb?query=duf285
mydata <- fread(file='DUF285_UniProtKB.tsv')


mydata$domain <- gsub(mydata$`Taxonomic lineage`, pattern = 'cellular organisms \\(no rank\\), ', replacement = '')
mydata$domain <- gsub(mydata$domain, pattern = ' .*', replacement = '')

mydata$phylum <- gsub(mydata$`Taxonomic lineage`, pattern = 'cellular organisms \\(no rank\\), ', replacement = '')
mydata$phylum <- sub(mydata$phylum, pattern = '\\(phylum\\), .*', replacement = '')
mydata[grep(mydata$phylum, pattern = 'Asgard'),]$phylum <- 'Asgard'
mydata[grep(mydata$phylum, pattern = 'DPANN'),]$phylum <- 'DPANN'
mydata[grep(mydata$phylum, pattern = 'TACK'),]$phylum <- 'TACK'
mydata$phylum <- sub(mydata$phylum, pattern = '\\(class\\), .*', replacement = '')
mydata$phylum <- gsub(mydata$phylum, pattern = '^.*), ', replacement = '')

#mydata <- mydata[mydata$domain!='Viruses',]
mydata <- mydata[mydata$domain!='Bacteria| environmental samples',]

############

#P. tricornutum has proteomes for 
mydata$Proteomes <- sub(mydata$Proteomes, pattern = 'UP000836788.*', replacement = 'UP000836788')
mydata$Proteomes <- sub(mydata$Proteomes, pattern = 'UP000000759.*', replacement = 'UP000000759')

mydata_agg <- aggregate(data=mydata, FUN=length, Entry ~ Proteomes + Organism + domain + phylum)
mydata_agg$Entry <- as.numeric(as.character(mydata_agg$Entry))
mydata_agg <- aggregate(data=mydata_agg, FUN=mean, Entry ~ Organism + domain + phylum)

mydata_agg <- mydata_agg[grep(x=mydata_agg$phylum, pattern = '\\(no rank\\)', invert = T),]

mydata_agg$phylum <- factor(mydata_agg$phylum, levels = 
gsub(x=sort(unique(paste(mydata_agg$domain, mydata_agg$phylum, sep='|'))), pattern = '^.*\\|', replacement = '')
)

mydata_agg$phylum <- factor(mydata_agg$phylum, levels = 
c("Asgard","Candidatus Thermoplasmatota ","DPANN","Euryarchaeota ","TACK","Actinomycetota ","Armatimonadota ","Bacillota ","Bacteroidota ","Balneolota ","Bdellovibrionota ","Caldisericota ","Calditrichota ","Campylobacterota ","candidate division CPR1 (phylum)","Candidatus Absconditabacteria ","Candidatus Gracilibacteria ","Candidatus Kaiserbacteria ","Candidatus Magasanikbacteria ","Candidatus Magasanikbacteria (phylum)","Candidatus Marinimicrobia ","Candidatus Moranbacteria ","Candidatus Muirbacteria ","Candidatus Parcubacteria ","Candidatus Peregrinibacteria ","Candidatus Poribacteria ","Candidatus Riflebacteria ","Candidatus Saccharibacteria ","Chlorobiota ","Chloroflexota ","Cyanobacteriota ","Deferribacterota ","Deinococcota ","Deltaproteobacteria ","Fibrobacterota ","Fusobacteriota ","Gemmatimonadota ","Ignavibacteriota ","Lentisphaerota ","Mycoplasmatota ","Myxococcota ","Planctomycetota ","Pseudomonadota ","Rhodothermota ","Spirochaetota ","Synergistota ","Thermodesulfobacteriota ","Verrucomicrobiota ","Annelida ","Arthropoda ", "Chordata ", "Basidiomycota ","Fornicata ","Mucoromycota ","Nematoda ","Oomycota ","Palpitomonas (genus)","Bigyra ","Bacillariophyta ", "Bolidophyceae ","Dictyochophyceae ","Pelagophyceae ","Pinguiophyceae ","Dinophyceae ","Haptophyta ","Chlorophyta ","Prasinodermophyta ","Rhodophyta ","Streptophyta ","Nucleocytoviricota ","Uroviricota ")
)

# Calculate sample sizes per phylum
library(tidyverse)
sample_sizes <- mydata_agg %>%
  group_by(phylum) %>%
  summarise(n = n(), .groups = 'drop')

# Add a new column for the label that includes 'n='
sample_sizes$label <- paste0("n=", sample_sizes$n)

# Calculate the maximum y position after log2 transformation for y position
max_y_position <- max(log2(mydata_agg$Entry), na.rm = TRUE) + 1  # Adjust the offset as needed

# Create a new data frame for the labels to ensure they're placed at the end of each plot
annotation_data <- sample_sizes %>%
  mutate(y_position = max_y_position)

# Now create the plot with the labels
figS17 <- ggplot(mydata_agg, aes(x = phylum, y = Entry, fill = domain)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal() +
  scale_y_continuous(trans = 'log2') +
  labs(title = "DUF285 (PF03382)",
       x = "Phylum",
       y = "# proteins from a proteome matched by at least one entry from DUF285 (log2 scale)",
       fill = 'Domain') +
  coord_flip() +
  geom_text(data = annotation_data, aes(x = phylum, y = y_position, label = label), 
            vjust = 0.5, hjust = 0, size = 3, color = "red", inherit.aes = FALSE)

figS17
