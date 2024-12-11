

# Guillaume Bernard
# guillaum3.bernard@gmail.com

library("WGCNA")
allowWGCNAThreads() 
options(stringsAsFactors = FALSE);

args = commandArgs(trailingOnly=TRUE)
getwd();
workingDir = "/Users/6ui11aum3/Documents/ibens/juan_v9_diatoms/V4/wgcna_diatoms_V4_pooledsize_filtered"
setwd(workingDir); 
lnames = load(file = "diatoms_V4_pooledsize_dataInput.RData")

#great at 12
softP = 12
mergeCutH = 0.75
detectCutH = 0.999

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, RsquaredCut = 0.9,verbose = 5, networkType = "signed")

# Plot the results:
sizeGrWindow(9, 5)
pdf("samples_soft_treshold_fit.pdf")
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()

pdf("samples_soft_treshold_connectivity.pdf")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
datExpr0[] <- lapply(datExpr0, as.numeric)
net = blockwiseModules(datExpr0, power = softP,
                       TOMType = "signed", minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = mergeCutH,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
dendro=paste("cluster_dendogram_pow",softP,"_merge",mergeCutH,"minMod20.pdf",sep="")
pdf(dendro)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
end_data=paste("diatoms_pfams_networkConstruction-pow",softP,"_merge",mergeCutH,"minMod20.RData",sep="")
save(MEs, moduleLabels, moduleColors, geneTree, softP, detectCutH, mergeCutH, traitColors, datTraits,
     file = end_data)
