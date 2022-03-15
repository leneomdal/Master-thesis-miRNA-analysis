# CREATE HEATMAPS USING ONLY TOP MIRNA FROM NEW VOOM ANALYSIS
source("Code//voom-analysis.R")
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(bioDist)
library(edgeR)
library(stats)
library(gplots)
library(dendextend)

#create df for annotations to use in heatmap
groups_sex = ifelse(metadata.df$sex == 1, "boy", "girl")
probiotic.groups = ifelse(metadata.df$probiotic == 1, "Probiotic", "Placebo")
ad.groups = ifelse(metadata.df$ad == 1, "AD", "No AD")
mad.groups = ifelse(metadata.df$matatopy == 1, "Maternal AD", "No maternal AD")
sample.col = data.frame(Supplement = probiotic.groups, Atopic.dermatitis = ad.groups, 
                        Maternal.atopy = mad.groups)
rownames(sample.col) = colnames(log.cpm)


# Plot 1-abs(cor) distance of top mirna from voom analysis
log.cpm.top.mirna = log.cpm[top.mirna,]

# drop one sample, sample "149" which is at place 39 !QUICK FIX!
log.cpm.top.mirna = log.cpm.top.mirna[,-39]




# calculate correlation between miRNAs and samples
cor.mirna.top = cor(t(log.cpm.top.mirna))
cor.samples.top = cor(log.cpm.top.mirna)

for(i in 1:ncol(log.cpm.top.mirna)){
  if(var(log.cpm.top.mirna[,i]) == 0){
    print(colnames(log.cpm.top.mirna)[i])
  }
}

probiotic.groups = ifelse(metadata.df$probiotic == 1, "Probiotic", "Placebo")
ad.groups = ifelse(metadata.df$ad == 1, "AD", "No AD")
mad.groups = ifelse(metadata.df$matatopy == 1, "Maternal atopy", "No maternal atopy")

sample.col = data.frame(Supplement = probiotic.groups, Atopic.dermatitis = ad.groups, 
                        Maternal.atopy = mad.groups)

pheatmap(log.cpm.top.mirna, scale = "row", border_color = NA, 
         clustering_distance_rows = as.dist(1-(cor.mirna.top)), 
         clustering_distance_cols = as.dist(1- (cor.samples.top)), 
         clustering_method = "ward.D2", annotation_col = sample.col, show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Heatmap of top miRNA using 1 - correlation distance",
         fontsize = 7)
