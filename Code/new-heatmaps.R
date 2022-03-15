source("Code//normalize-and-filter.R")
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(bioDist)
library(edgeR)
library(stats)
library(gplots)
library(dendextend)

# set colors
colors <- colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(255)

mirna.clust = as.dendrogram(hclust(dist(log.cpm)))
sample.clust = as.dendrogram(hclust(dist(t(log.cpm))))
colSide <- brewer.pal(9, "Set1")[groups_ad]
hclust.mirna = hclust(as.dist(1-abs(cor(t(log.cpm)))), method = "ward.D2")
as.dendrogram(hclust.mirna) %>% plot(horiz = TRUE)


groups_sex = ifelse(metadata.df$sex == 1, "boy", "girl")

probiotic.groups = ifelse(metadata.df$probiotic == 1, "Probiotic", "Placebo")
ad.groups = ifelse(metadata.df$ad == 1, "AD", "No AD")
mad.groups = ifelse(metadata.df$matatopy == 1, "Maternal atopy", "No maternal atopy")

sample.col = data.frame(Supplement = probiotic.groups, Atopic.dermatitis = ad.groups, 
                        Maternal.atopy = mad.groups)

rownames(sample.col) = colnames(log.cpm)


# Heatmap using euclidian distance
clust.m.euclid = hclust(dist(t(scale(t(log.cpm)))), method = "ward.D2")
clust.s.euclid = hclust(dist(scale(t(log.cpm))), method = "ward.D2")
View(t(scale(t(log.cpm))))


pheatmap(log.cpm,
         col = colors, scale = "row",
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D2", annotation_col = sample.col, show_rownames = FALSE,
         main = "Heatmap of log cpm values using euclidian distance", show_colnames = FALSE,
        fontsize = 7, border_color = NA)


# calculate correlation between miRNAs and samples
cor.mirna = cor(scale(t(log.cpm)))
cor.samples = cor(t(scale(t(log.cpm))))


# Heatmap using cor based disance !!!! 1-cor !!!!
clust.m.cor = hclust(as.dist(1-cor.mirna), method = "ward.D2")
clust.s.cor = hclust(as.dist(1-cor.samples), method = "ward.D2")

pheatmap(log.cpm, color = colors, scale = "row",  border_color = NA, 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         clustering_method = "ward.D2", annotation_col = sample.col, show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of log cpm values using 1 - correlation distance", fontsize = 7)


clust.m.abs.cor = hclust(as.dist(1-abs(cor.mirna)), method = "ward.D2")
clust.s.abs.cor = hclust(as.dist(1-abs(cor.samples)), method = "ward.D2")
#plot 1-abs(cor) distance
pheatmap(log.cpm, scale = "row", border_color = NA, 
         clustering_distance_rows = as.dist(1-abs(cor.mirna)), 
         clustering_distance_cols = as.dist(1- abs(cor.samples)), 
         clustering_method = "ward.D2", annotation_col = sample.col, show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of log cpm values using 1 - abs(correlation) distance",
         fontsize = 7)



