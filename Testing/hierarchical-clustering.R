# From end to end RNA-seq analysis with deseq

source("Code//load_count_data.R")
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(bioDist)
library(edgeR)
library(stats)

#Exploration of mean-variance relationship
# plot mean vs variance across samples
means = rowMeans(assays(dds)$counts)
sds = apply((assays(dds)$counts), 1, sd)
plot(means, sds^2)
# Plot first sample against second
plot(as.numeric(assays(dds)$counts[,1]), as.numeric(assays(dds)$counts[,2]))
# Plot mean variance after log transform
means.rlog = rowMeans(assay(rld))
sds.rlog = apply(assay(rld), 1, sd)
plot(means.rlog, sds.rlog)
# Plot mean variance after log cpm transform
mean.log.cpm = rowMeans(log.cpm)
sds.log.cpm = apply(log.cpm, 1, sd)
plot(mean.log.cpm, sds.log.cpm)

# Heatmap 
mirna.clust = as.dendrogram(hclust(dist(log.cpm)))
sample.clust = as.dendrogram(hclust(dist(t(log.cpm))))
colSide <- brewer.pal(9, "Set1")[groups_p]
heatmap(log.cpm, scale = NULL, distfun = dist, ColSideColors = colSide)


# set colors
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


# Intial heatmaps



# First using rlog values

# Calculate euclidian sample distances
euclid.dist.samples.rld = stats::dist(t(assay(rld)))
euclid.dist.samples.rld.mat = as.matrix(euclid.dist.samples.rld)
rownames(euclid.dist.samples.rld.mat) = metadata.df$probiotic
colnames(euclid.dist.samples.rld.mat) <- NULL
pheatmap(euclid.dist.samples.rld.mat,
         clustering_distance_rows=euclid.dist.samples.rld,
         clustering_distance_cols=euclid.dist.samples.rld, col=colors)



# Correlation distance, samples
cor.dist.samples.rld <- cor.dist(t(assay(rld)))
cor.dist.samples.rld.matrix <- as.matrix(cor.dist.samples.rld)
rownames(cor.dist.samples.rld.matrix) = NULL
colnames(cor.dist.samples.rld.matrix) = groups_ad
pheatmap(cor.dist.samples.rld.matrix,
         clustering_distance_rows=cor.dist.samples.rld,
         clustering_distance_cols=cor.dist.samples.rld,
         col=colors)

# Euclidian distance, miRNAs
euclid.dist.mirna.rld <- dist(assay(rld))
euclid.dist.mirna.rld.mirna <- as.matrix(euclid.dist.mirna.rld)
rownames(euclid.dist.mirna.rld.mirna) = NULL
colnames(euclid.dist.mirna.rld.mirna) = rownames(rld)
pheatmap(euclid.dist.mirna.rld.mirna,
         clustering_distance_rows=euclid.dist.mirna.rld,
         clustering_distance_cols=euclid.dist.mirna.rld,
         col=colors)



# Correlation distance, miRNAs
#cor.dist.mirna.rld <- cor.dist(counts(dds), abs = TRUE)
cor.dist.mirna.rld <- cor.dist(assay(rld))
cor.dist.mirna.rld.matrix <- as.matrix( cor.dist.mirna.rld )
rownames(cor.dist.mirna.rld.matrix) <- NULL
colnames(cor.dist.mirna.rld.matrix) <- rownames(rld)
pheatmap(cor.dist.mirna.rld.matrix,
         clustering_distance_rows = cor.dist.mirna.rld,
         clustering_distance_cols = cor.dist.mirna.rld,
         col=colors, scale = "none")



# USING LOG CPM VALUES


# Euclidian distance, samples
euclid.dist.samples.lcpm = dist(t(log.cpm), method = "euclidian")
euclid.dist.samples.lcpm.matrix = as.matrix(euclid.dist.samples.lcpm)
colnames(euclid.dist.samples.lcpm.matrix) = groups_ad
pheatmap(euclid.dist.samples.lcpm.matrix,
         clustering_distance_row = euclid.dist.samples.lcpm,
         clustering_distance_cols = euclid.dist.samples.lcpm,
         col=colors, clustering_method = "complete")

# Correlation distance, samples
cor.dist.samples.lcpm = cor.dist(t(log.cpm), abs = TRUE)
cor.dist.samples.lcpm.mat = as.matrix(cor.dist.samples.lcpm)
colnames(cor.dist.samples.lcpm.mat) = groups_p
pheatmap(cor.dist.samples.lcpm.mat,
         clustering_distance_rows = cor.dist.samples.lcpm,
         clustering_distance_cols = cor.dist.samples.lcpm,
         col=colors, scale = "none",  clustering_method = "complete")


# Euclidian distance, miRNAs
euclid.dist.mirna.lcpm = dist(log.cpm, method = "euclidian")
euclid.dist.mirna.lcpm.matrix = as.matrix(euclid.dist.mirna.lcpm)
colnames(euclid.dist.mirna.lcpm.matrix) = rownames(log.cpm)
rownames(euclid.dist.mirna.lcpm.matrix)  = NULL
pheatmap(euclid.dist.mirna.lcpm.matrix,
         clustering_distance_row = euclid.dist.mirna.lcpm,
         clustering_distance_cols = euclid.dist.mirna.lcpm,
         col=colors, clustering_method = "complete")

# Correlation distance, miRNAs
cor.dist.mirna.lcpm = cor.dist(log.cpm, abs = TRUE)
cor.dist.mirna.lcpm.mat = as.matrix(cor.dist.mirna.lcpm)
colnames(cor.dist.mirna.lcpm.mat) = rownames(log.cpm)
pheatmap(cor.dist.mirna.lcpm.mat,
         clustering_distance_rows = cor.dist.mirna.lcpm,
         clustering_distance_cols = cor.dist.mirna.lcpm,
         col=colors, scale = "none",  clustering_method = "complete")




# Poisson distance sample clustering
# The PoissonDistance function takes the original count matrix (not normalized, as it normalizes counts internally) with samples as rows instead of columns, so we need to transpose the counts in dds.
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- rld$probiotic
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)



#MODEL BASED CLUSTERING

library(MBCluster.Seq)
countdata.MB = RNASeq.Data(Count = count.df, Normalize = log(scalar), 
                           Treatment = metadata.df$probiotic, GeneID = rownames(count.df))








