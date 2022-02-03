
source("Code//load_count_data.R")
library(bioDist)
library(edgeR)
#Filter and normalize


#Make dge object
dge <- DGEList(counts=count.df)
#Convert to cpm
cpm <- cpm(dge)
# Filter out only those with all 0 counts
keep = rowSums(dge$counts==0)!=60
count.df = count.df[keep, ]
dim(count.df)

#Hierarchical clustering

# Modify the location based on your filesystem
load("Testing/pca-examples.Rdata")

# We will work with nyt.frame
nyt_data = nyt.frame
hc.complete = hclust(dist(nyt_data[,-1]), method = "complete")
plot(hc.complete)



# cluster considering samples
hc.samples.abs = hclust(as.dist(1- abs(cor(t(as.matrix(scale(t(count.df))))))))
hc.samples = hclust(as.dist(1- cor(t(as.matrix(scale(t(count.df)))))))
hc.samples = hclust(cor.dist(as.matrix(scale(t(count.df))), abs = FALSE))
plot(hc.samples.abs)

#clusters considering miRNAs
hc.miRNA = hclust(as.dist(1 - cor(as.matrix(t(scale(count.df))))))
plot(hc.miRNA)


# Add groups to count matrix
#Make groups
treatment = ifelse(metadata.df$probiotic ==1, "P", "nP")
outcome = ifelse(metadata.df$ad == 1, "AD", "nAD")
groups_p = as.factor(treatment)
groups = as.factor(paste0(treatment,"_", outcome))
groups_ad = as.factor(outcome)
count.df.g = count.df
count.df.g["AD",] = groups_ad
count.df.g["P", ] = groups_p
count.df.g["groups", ] = groups
dim(count.df.g)




#PCA
#eksempel
nyt_pca = prcomp(nyt_data[,-1])
biplot(nyt_pca, scale=0)
nyt_loading = nyt_pca$rotation[, 1:2]

#bør det ikke være absoluttverdi når man sorterer her?
informative_loadings = rbind(
  head(nyt_loading[order(nyt_loading[,1], decreasing = TRUE),]),
  head(nyt_loading[order(nyt_loading[,2], decreasing = TRUE),])
)
biplot(x = nyt_pca$x[,1:2], y= informative_loadings, scale=0)



#miRNA data
counts.scaled = scale(t(count.df))
expr_pca = prcomp(counts.scaled)
expr_loading = expr_pca$rotation[,1:2]
informative_loadings = rbind(
  head(expr_loading[order(abs(expr_loading[,1]), decreasing = TRUE),]),
  head(expr_loading[order(abs(expr_loading[,2]), decreasing = TRUE),])
)
biplot(x = expr_pca$x[,1:2], y= informative_loadings, scale=0)
