
library(ggplot2)
library(cowplot)


################### distribution of top 100 diatom Pfams among size fractions and ocean regions
source("organizer_fig9/top100_metaT_diatoms.R")
p
rm( "df" , "p", "OR_abundance", "pfam_abundance",  "size_abundance" , "df_ranked")

################# taxonomic distribution among eukaryotic phytoplankton of the top 100 diatom Pfams
source("organizer_fig9/top100_metaT_diatoms_vs_otherPhyto.R")
rm('final_data', 'pfam_data', 'pfam_data_merged', 'top100Pfams', 'top100', 'mycolors')

Fig9=plot_grid(p_abundance, p_size, pOR, p_phyto, ncol=4, rel_widths = c(3,1,1,1), align="h")
Fig9



