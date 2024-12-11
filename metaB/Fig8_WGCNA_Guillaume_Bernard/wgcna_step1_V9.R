

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
workingDir = "/Users/6ui11aum3/Documents/ibens/juan_v9_diatoms/V9"
setwd(workingDir); 
#file_mat = args[2]
file_mat = "V9_abundance_matrix.pooledsizes.csv"

#res_dir = args[3];
res_dir = "wgcna_diatoms_V9_filtered_pooled"
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
traitData = read.csv("V9_environment.pooledsizes_v2.csv");
dim(traitData)
names(traitData)

# necessary transformations (example: latitude to abs latitude)
traitData <- traitData %>% rename(Latitude = lat)
traitData <- traitData %>% rename(Longitude = long)
traitData$Latitude <- abs(traitData$Latitude)
traitData <- traitData %>% rename(Abs_Latitude = Latitude)


#traitData$depth <- replace(traitData$depth, traitData$depth=="SRF", 1)
#traitData$depth <- replace(traitData$depth, traitData$depth=="MIX", 2)
#traitData$depth <- replace(traitData$depth, traitData$depth=="DCM", 2)

# remove columns that hold information we do not need.
allTraits = traitData[, c(1:51)];

dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

pfamSamples = rownames(datExpr0);
traitRows = match(pfamSamples, allTraits$sample);
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
save(datExpr, datExpr0, traitColors, datTraits, file = paste(res_dir,"diatoms_V9_pooled_dataInput.RData",sep="/"))


###### BONUS ##### spearman correlation of env matrix #####

env_matrix <- datTraits %>% select(-Longitude)
env_matrix <- env_matrix %>% select(-NPP.uqartakuvik)

mat_1 <-round(cor(env_matrix, method = "spearman", use = 'pairwise.complete.obs'),2)
mat_1

# Exp 1 #
library("GGally")

ggcorr(env_matrix, method = c("pairwise", "spearman"),
       nbreaks = NULL, digits = 2, low = "steelblue", mid = "white", high = "darkred",
       geom = "tile", label = FALSE,
       label_size = 3, label_round = 2, label_alpha = FALSE,
       hjust = 0.75, size = 3, color = "grey50", layout.exp = 5)

# Exp 2 #

library(gapmap)
library(RColorBrewer)
RdBu = rev(brewer.pal(11, name="RdBu"))
RdYlBu = rev(brewer.pal(11, name="RdYlBu"))


env_matrix2 <- env_matrix %>% select(-PAR.PC)
env_matrix2 <- env_matrix2 %>% select(-acCDOM.atsushi)
env_matrix2 <- env_matrix2 %>% select(-sCDOM.atsushi)
env_matrix2 <- env_matrix2 %>% select(-Ice.free.period)


row_dist <- as.dist(1-cor(env_matrix, method = "spearman", use = 'pairwise.complete.obs'))
col_dist <- as.dist(1-cor(env_matrix, method = "spearman", use = 'pairwise.complete.obs'))

col_hc <- hclust(col_dist, method = "complete")
row_hc <- hclust(row_dist, method = "complete")
col_d <- as.dendrogram(col_hc)
row_d <- as.dendrogram(row_hc)

gapmap(m = as.matrix(env_matrix2), d_row = rev(row_d), d_col = col_d, ratio = 0, verbose=FALSE, col=RdBu,
       label_size=2, v_ratio= c(0.1,0.8,0.1), h_ratio=c(0.1,0.8,0.1))

row_data <- gap_data(d= row_d, mode = "quantitative", mapping="exponential", ratio=0.3, scale= 0.5)
dend <- gap_dendrogram(data = row_data, leaf_labels = TRUE, rotate_label = TRUE)
dend + theme(axis.ticks.length= grid::unit(0,"lines") )+ theme(axis.ticks.margin = grid::unit(-0.8, "lines"))


# Exp 3 #

library(gplots)

mycolors <- colorRampPalette(rev(brewer.pal(11,'RdYlBu')))
heatmap.2(mat_1, trace = "none", col = mycolors)



