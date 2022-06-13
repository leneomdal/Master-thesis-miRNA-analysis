# CREATE HEATMAPS USING ONLY TOP MIRNA FROM VOOM ANALYSIS
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
#log.cpm.top.mirna = log.cpm.top.mirna[,-39]


# calculate correlation between miRNAs and samples
cor.mirna.top = cor(scale(t(log.cpm.top.mirna)), method = "pearson")
cor.samples.top = cor(t(scale(t(log.cpm.top.mirna))), method = "pearson")

# Without scaling
#cor.mirna.top = cor(t(log.cpm.top.mirna))
#cor.samples.top = cor(log.cpm.top.mirna)

# Only centering
#cor.mirna.top = cor(scale(t(log.cpm.top.mirna), scale = FALSE))
#cor.samples.top = cor(t(scale(t(log.cpm.top.mirna), scale = FALSE)))



# data frame for annotations
probiotic.groups = ifelse(metadata.df$probiotic == 1, "Probiotic", "Placebo")
ad.groups = ifelse(metadata.df$ad == 1, "AD", "No.AD")
mad.groups = ifelse(metadata.df$matatopy == 1, "Maternal.atopy", "No.maternal.atopy")
sex.groups = ifelse(metadata.df$sex == 1, "boy", "girl")
sample.col = data.frame(Supplement = as.factor(probiotic.groups), Atopic.dermatitis = as.factor(ad.groups), 
                        Maternal.atopy = as.factor(mad.groups))
rownames(sample.col) = colnames(log.cpm)

#Set colors for annotations
ann.colors = list(
  Supplement = c(Placebo = "#BDFCC9", Probiotic = "#03A89E"),
  Atopic.dermatitis = c(AD = "#8B4789",No.AD = "#87CEEB"),
  Maternal.atopy = c(Maternal.atopy = "#FFA07A",No.maternal.atopy = "#1B9E77")
)


# 1-abs(corr)
pheatmap(log.cpm.top.mirna, scale = "row", border_color = NA, 
         clustering_distance_rows = as.dist(1-abs(cor.mirna.top)), 
         clustering_distance_cols = as.dist(1- abs(cor.samples.top)), 
         clustering_method = "ward.D2", annotation_col = sample.col, annotation_colors = ann.colors, show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Heatmap of top miRNA using 1 - abs(correlation) distance",
         fontsize = 7)
# 1-corr
pheatmap(log.cpm.top.mirna, scale = "row", border_color = NA, 
         clustering_distance_rows = as.dist(1-(cor.mirna.top)), 
         clustering_distance_cols = as.dist(1- (cor.samples.top)), 
         clustering_method = "ward.D2", annotation_col = sample.col, annotation_colors = ann.colors, show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Heatmap of top miRNA using 1 - correlation distance",
         fontsize = 7)
# Euclidean
pheatmap(log.cpm.top.mirna, scale = "row", border_color = NA, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2", annotation_col = sample.col, annotation_colors = ann.colors, show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Heatmap of top miRNA using euclidean distance",
         fontsize = 7)


# Create hclust object (current: correlation based distance)
clust.s.top = hclust(as.dist(1-(cor.samples.top)), method = "ward.D2")
plot(clust.s.top)
