
# Guillaume Bernard
# guillaum3.bernard@gmail.com

library(WGCNA);
library(dplyr);
library(readr)
library(janitor)

args = commandArgs(trailingOnly=TRUE)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = args[1];
workingDir = "/Users/6ui11aum3/Documents/ibens/juan_v9_diatoms/V4"
setwd(workingDir); 
#file_mat = args[2]
file_mat = "V4_abundance_matrix.pooledsizes.csv"

#res_dir = args[3];
res_dir = "wgcna_diatoms_V4_pooledsize_filtered"
dir.create(res_dir);
##### LOAD DATA

data = read.csv(file_mat)
datExpr0 <- as.data.frame(t(data[,-c(1)]));
     
names(datExpr0) = data$ASV;
rownames(datExpr0) = names(data)[-c(1)];
#rownames(datExpr0) <- lapply(rownames(datExpr0) , function(x){substring(x, 2)})

#################
#filter genes that appears in more than 2 stations
my_filter = function(col) {
  sum(col==0)<length(col)-2
}

#filter genes that appears more than 3 times
my_filter2 = function(col) {
  sum(col)>3
}

datExpr1=datExpr0[,apply(datExpr0,2,my_filter)]
datExpr0=datExpr1[,apply(datExpr1,2,my_filter2)]
#################

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
#pdf(paste(res_dir,"sampleClustering_pgcna.pdf",sep="/"))
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#dev.off()
#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

# Plot a line to show the cut
#abline(h = 20000, col = "red");
# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr0[keepSamples, ]
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#datExpr <- log10(datExpr)

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

#traitData = read.csv(args[4]);
traitData = read.csv("V4_environment.pooledsizes_v2.csv");
dim(traitData)
names(traitData)

# necessary transformations (example: latitude to abs latitude)
traitData$Latitude <- abs(traitData$Latitude)
traitData <- traitData %>% rename(Abs_Latitude = Latitude)


#traitData$Depth <- replace(traitData$Depth, traitData$Depth=="SRF", 1)
#traitData$Depth <- replace(traitData$Depth, traitData$Depth=="MIX", 2)
#traitData$Depth <- replace(traitData$Depth, traitData$Depth=="DCM", 2)
traitData$Depth.nominal <- replace(traitData$Depth.nominal, traitData$Depth.nominal=="0-100", 50)
traitData$Depth.nominal <- replace(traitData$Depth.nominal, traitData$Depth.nominal=="10-100", 50)

# remove columns that hold information we do not need.
allTraits = traitData[, c(1:58)];

dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

pfamSamples = rownames(datExpr0);
traitRows = match(pfamSamples, allTraits$Sample.mat);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
datTraits <- datTraits %>% mutate_if(is.character, as.numeric)
traitColors = numbers2colors(datTraits, signed = TRUE);

# Plot the sample dendrogram and the colors underneath.
pdf(paste(res_dir,"samples_trait_heatmap.pdf",sep="/"))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

#save(datExpr, datTraits, file = paste(res_dir,args[5],sep="/"))
save(datExpr, datExpr0, traitColors, datTraits, file = paste(res_dir,"diatoms_V4_pooledsize_dataInput.RData",sep="/"))


# Exp 3 #

datTraits <- datTraits %>% mutate_if(is.character, as.numeric)

env_matrix <- datTraits %>% select(-Longitude)
env_matrix <- env_matrix %>% select(-NPP.uqartakuvik)
env_matrix <- env_matrix %>% select(-Depth.euph.zone)
env_matrix <- env_matrix %>% select(-Depth.Mixed.Layer)
env_matrix <- env_matrix %>% select(-Depth.chloro.max)
env_matrix <- env_matrix %>% select(-Depth.max.Brunt.Väisälä)
env_matrix <- env_matrix %>% select(-Depth.Max.O2)
env_matrix <- env_matrix %>% select(-Depth.Min.O2)
env_matrix <- env_matrix %>% select(-Depth.nitracline)

mat_1 <-round(cor(env_matrix, method = "spearman", use = 'pairwise.complete.obs'),2)
mat_1

library(RColorBrewer)
RdBu = rev(brewer.pal(11, name="RdBu"))
RdYlBu = rev(brewer.pal(11, name="RdYlBu"))

library(gplots)

mycolors <- colorRampPalette(rev(brewer.pal(11,'RdYlBu')))
heatmap.2(mat_1, trace = "none", col = mycolors)

