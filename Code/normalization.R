rm(list=ls()) 
library(edgeR)
library(xtable)
library(RColorBrewer)
library(stats)
source("Code//load_count_data.R")


#Make dge object
dge <- DGEList(counts=count.df)

#Make groups
treatment = ifelse(metadata.df$probiotic == 1, "P", "nP")
outcome = ifelse(metadata.df$ad == 1, "AD", "nAD")

groups_p = as.factor(treatment)
groups = as.factor(paste0(treatment,"_", outcome))
groups_ad = as.factor(outcome)
groups_mad = as.factor(ifelse(metadata.df$matatopy == 1, "mAD", "no mAD"))



#Convert to cpm and log cpm 
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)
summary(lcpm)


# TABLE: not expressed miRNAs
table.zero.expression = table(rowSums(dge$counts==0)==60)

# Filtering out miRNAs with cpm value less than 500 in 10% or less of the samples
#keep = rowSums(cpm(dge)>500)>=6
#dge$counts = dge$counts[keep,]
#dim(dge)


#Second way of filtering out low expressed miRNA
#keep.exprs <- filterByExpr(dge, group = groups)
#dge.g <- dge[keep.exprs,, keep.lib.sizes=FALSE]
#dim(dge.g$counts)

# Filter out only those with all 0 counts
keep = rowSums(dge$counts==0)!=60
sum(keep)
dge$counts = dge$counts[keep,]
dim(dge$counts)
dim(count.df)


#Optional normalization. Not used here  https://support.bioconductor.org/p/69433/
dge = calcNormFactors(dge, method = "TMM")

count.matrix = as.matrix(dge$counts)
View(count.matrix)
count.scaled = scale(t(count.matrix))


View(cor(t(count.matrix)))
heatmap(count.matrix, Rowv = 'dendrogram')
?heatmap

cU = cor(t(count.scaled))

heatmap(cU, Rowv = FALSE, symm = TRUE, distfun = function(c) as.dist(1 - c), keep.dendro = TRUE, scale = NULL)
View(metadata.df)

print(metadata.df$samples[metadata.df$ad == 1])

plot(hclust(cor(t(count.scaled))))