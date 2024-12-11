#!/usr/bin/env Rscript --vanilla
#The --vanilla on the end, tells Rscript to run without saving or restoring anything

rm(list=ls())

args = commandArgs(trailingOnly=TRUE)

library(tidyr)
library(tidyverse)
library(data.table)


calculador <- function(MYtaxon, Mymeta){

pfam <- fread(file=
 paste(args[1], ".matou-v1.5.Pfam", sep=""),
 sep="\t", header=T)

pfam <- pfam %>% select(geneID, pfamAcc)



meta_long <- fread(file=
 paste(args[1], ".matou-v1.5.", Mymeta, ".occurrences.tsv", sep="")
)


colnames(meta_long) <- c("geneID", "sample", "rpkm")

meta_long_pfam <- merge(x=pfam, y=meta_long, all.y=T, by="geneID", allow.cartesian=TRUE)
rm(meta_long)

meta_long_pfam_agg <- aggregate(data=meta_long_pfam, FUN=sum, rpkm ~ pfamAcc + sample)
rm(meta_long_pfam)

meta_long_pfam_agg$taxon <- args[1]

meta_long_pfam_agg$meta <- Mymeta

write.table(meta_long_pfam_agg, file=paste(args[1], ".matou-v1.5.Pfam.", Mymeta, ".tsv", sep=""), sep="\t", row.names=F, quote=F)
rm(meta_long_pfam_agg)

}



calculador(args[1], "metaG")
calculador(args[1], "metaT")

