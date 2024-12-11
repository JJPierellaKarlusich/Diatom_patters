

# Guillaume Bernard
# guillaum3.bernard@gmail.com

library("WGCNA")
library("janitor")
library(dplyr)
library(ape)

args = commandArgs(trailingOnly=TRUE)
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = args[1];
workingDir = "/Users/6ui11aum3/Documents/ibens/juan_v9_diatoms/V4/wgcna_diatoms_V4_pooledsize_filtered"
setwd(workingDir); 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
allowWGCNAThreads() 

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

# Load network data saved in the second part.
#lnames = load(file = args[3]);
lnames = load(file = "diatoms_pfams_networkConstruction-pow12_merge0.75minMod20.RData")
lnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#Change latitude to abs values
#datTraits$latitude <- abs(datTraits$latitude)

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================



plot_mod = paste("Module-trait_relationships_pow",softP,"_merge",mergeCutH,"minMod20.pdf",sep="")
sizeGrWindow(40,24)
pdf(plot_mod)

# Will display correlations and their p-values

moduleTraitPvalue2 <- matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method='fdr'),ncol=57)
rownames(moduleTraitPvalue2) <- rownames(moduleTraitPvalue)
colnames(moduleTraitPvalue2) <- colnames(moduleTraitPvalue)

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue2, 1), ")", sep = "");

write.csv(signif(moduleTraitCor, 2), "V4_module_envVariable_correlations_matrix.csv")
write.csv(signif(moduleTraitPvalue2, 1), "V4_module_envVariable_pvalues_matrix.csv")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 8, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.05,
               cex.lab.x = 0.3,
               cex.lab.y = 0.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

dir_output = "gene_MM_GS_info_pow12_merge0.25"
dir.create(dir_output)
setwd(dir_output); 
for (i in colnames(datTraits)){
  
  trait_to_analyse = i
  # Define variable weight containing the weight column of datTrait
  weight = as.data.frame(datTraits[[i]]);
  names(weight) = "weight"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  
  names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
  names(GSPvalue) = paste("p.GS.", names(weight), sep="");
  
  #=====================================================================================
  #
  #  Code chunk 5
  #
  #=====================================================================================
  for (j in modNames){
    module = j
    column = match(module, modNames);
    moduleGenes = moduleColors==module;
    pdf(paste(i,"_",j,"_MMvsGS.pdf",sep=""))
    sizeGrWindow(7, 7);
    par(mfrow = c(1,1));
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for ", i, sep=""),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    dev.off()
    #=====================================================================================
    #
    #  Code chunk 6
    #
    #=====================================================================================
    
    names(datExpr0)
    
    #=====================================================================================
    #
    #  Code chunk 7
    #
    #=====================================================================================
    
    names(datExpr)[moduleColors==j]
  }
  #annot = read.csv(file = "GeneAnnotation.csv");
  #dim(annot)
  #names(annot)
  probes = names(datExpr0)
  #probes2annot = match(probes, annot$substanceBXH)
  # The following is the number or probes without annotation:
  #sum(is.na(probes2annot))
  # Should return 0.
  
  #=====================================================================================
  #
  #  Code chunk 9
  #
  #=====================================================================================
  
  # Create the starting data frame
  geneInfo0 = data.frame(substanceBXH = probes,
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  # Order modules by their significance for weight
  modOrder = order(-abs(cor(MEs, weight, use = "p")));
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
  geneInfo = geneInfo0[geneOrder, ]
  
  #=====================================================================================
  #
  #  Code chunk 10
  #
  #=====================================================================================
  
  write.csv(geneInfo, file = paste("genesInfo_",i,".csv",sep = ""))
}


#####     module heatmap and the eigengene
datME = MEs0
signif(cor(datME, use="p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")

my_tree <- as.phylo(hclustdatME) 
write.tree(phy=my_tree, file="V4_module_clustering_tree.newick")

sizeGrWindow(8,9)
plotMEpairs(datME)

sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
for (j in modNames){

  which.module=j;
  plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),nrgcols=30,rlabels=T,
          clabels=T,rcols=which.module,
          title=which.module )
}

##### Relating modules and module eigengenes to external data
datME = MEs0
modNames = substring(names(MEs), 3)

count = 0;
sizeGrWindow(8,7);
for (j in modNames){

  which.module=j
  ME1=datME[, paste("ME",which.module, sep="")]
  #par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
  #plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
  #        nrgcols=30,rlabels=F,rcols=which.module,
  #        main=which.module, cex.main=2)
  #par(mar=c(5, 4.2, 0, 0.7))
  #barplot(ME1, col=which.module, main="", cex.main=2,
   #       ylab="eigengene expression",xlab="array sample", names.arg = rownames(datExpr), las=2, cex.names=0.3)
  if (count==0){
    combine_info <- data.frame( mod = ME1)
    names(combine_info)[names(combine_info) == "mod"] <- which.module
  }
  else{
    combine_info <- data.frame(combine_info,mod = ME1)
    names(combine_info)[names(combine_info) == "mod"] <- which.module
  }
  count = count + 1;
}

rownames(combine_info) = rownames(datExpr0)

combine_info <- cbind(combine_info, longitude = datTraits$Longitude)
combine_info <- cbind(combine_info, latitude = datTraits$Abs_Latitude)
combine_info <- cbind(combine_info, sample = rownames(datExpr0))

rm(pfam_info)
pfam_info <- datExpr0
pfam_info<- tibble::rownames_to_column(pfam_info, "Sample")
pfam_info <- pfam_info %>% adorn_totals("col", name = "Total count")
n <- ncol(pfam_all)
pfam_info2 <- pfam_info[,c(1,n)]

count = 0;
for (j in modNames){
  which.module=j
  pfam_info <- datExpr0[,moduleColors==which.module ]
  pfam_info<- tibble::rownames_to_column(pfam_info, "sample")
  pfam_module <- data.frame(OTU = colnames(datExpr), module = moduleColors)
  pfam_info <- pfam_info %>% adorn_totals("col", name = "M1count")
  if (count == 0){
    pfam_all <- pfam_info[,c("sample","M1count")]
    names(pfam_all)[names(pfam_all) == "M1count"] <- paste(which.module,"count",sep="")
  }
  else{
    pfam_all <- cbind(pfam_all, M1count = pfam_info[,c("M1count")])
    names(pfam_all)[names(pfam_all) == "M1count"] <- paste(which.module,"count",sep="")
  }
  count = count +1;

}

final_df <- cbind(combine_info[,c("longitude","latitude")],pfam_all[,c(1:n)])
final_df2 <- cbind(final_df, total = pfam_info2[,c("Total count")])

write.csv(pfam_module,"OTU_module_pow12_merge0.25_minMod40.csv",row.names = FALSE)
