source("Code//normalize-and-filter.R")


# K-means clustering on log cpm values

# scale by miRNA
scaled.log.cpm.values = scale(t(log.cpm)) 
k.m = kmeans(scaled.log.cpm.values, 2)
lcpm.km.cluster = k.m$cluster
?kmeans

# non-scaled values
log.cpm.kmeans = kmeans(t(log.cpm), 2)
log.cpm.km.cluster = log.cpm.kmeans$cluster

# K-means clustering on deseqs rlog transformed data

# scale by miRNA
scaled.rld.values = scale(t(assay(rld)))
kmeans.deseq.norm = kmeans(scaled.rld.values, 2)
rld.km.cluster = kmeans.deseq.norm$cluster

# not scaled values
rld.kmeans = kmeans(t(assay(rld)), 2)
rld.kmeans.cluster = rld.kmeans$cluster

# Compute PCA to view results

#counts.scaled = scale(t(count.df))
#log cpm
expr_pca = prcomp(t(log.cpm))
expr_loading = expr_pca$rotation[,1:2]
informative_loadings = rbind(
  head(expr_loading[order(abs(expr_loading[,1]), decreasing = TRUE),]),
  head(expr_loading[order(abs(expr_loading[,2]), decreasing = TRUE),])
)
biplot(x = expr_pca$x[,1:2], y= informative_loadings, scale=0)
plot(expr_pca$x[,1:2], type = "n")
points(expr_pca$x[groups_p == "P", 1:2], pch = "P", col = (log.cpm.km.cluster+1)[groups_p == "P"]) 
points(expr_pca$x[groups_p == "nP", 1:2], pch = "N", col = (log.cpm.km.cluster+1)[groups_p == "nP"])



# LCPM scaled
expr_pca = prcomp(scaled.log.cpm.values)
expr_loading = expr_pca$rotation[,1:2]
informative_loadings = rbind(
  head(expr_loading[order(abs(expr_loading[,1]), decreasing = TRUE),]),
  head(expr_loading[order(abs(expr_loading[,2]), decreasing = TRUE),])
)
biplot(x = expr_pca$x[,1:2], y= informative_loadings, scale=0)
plot(expr_pca$x[,1:2], type = "n")
points(expr_pca$x[groups_p == "P", 1:2], pch = "P", col = (lcpm.km.cluster+1)[groups_p == "P"]) 
points(expr_pca$x[groups_p == "nP", 1:2], pch = "N", col = (lcpm.km.cluster+1)[groups_p == "nP"])
expr_loading


#RLD

rld_pca = prcomp(scaled.rld.values)
rld_loading = rld_pca$rotation[,1:2]
rld_informative_loadings = rbind(
  head(rld_loading[order(abs(rld_loading[,1]), decreasing = TRUE),]),
  head(rld_loading[order(abs(rld_loading[,2]), decreasing = TRUE),])
)
biplot(x = rld_pca$x[,1:2], y= rld_informative_loadings, scale=0)
plot(rld_pca$x[,1:2], type = "n")
points(rld_pca$x[groups_p == "P", 1:2], pch = "P", col = (rld.km.cluster+1)[groups_p == "P"]) 
points(rld_pca$x[groups_p == "nP", 1:2], pch = "N", col = (rld.km.cluster+1)[groups_p == "nP"])

# Unscaled

rld_pca = prcomp(t(assay(rld)))
rld_loading = rld_pca$rotation[,1:2]
rld_informative_loadings = rbind(
  head(rld_loading[order(abs(rld_loading[,1]), decreasing = TRUE),]),
  head(rld_loading[order(abs(rld_loading[,2]), decreasing = TRUE),])
)
biplot(x = rld_pca$x[,1:2], y= rld_informative_loadings, scale=0)
plot(rld_pca$x[,1:2], type = "n")
points(rld_pca$x[groups_p == "P", 1:2], pch = "P", col = (rld.kmeans.cluster+1)[groups_p == "P"]) 
points(rld_pca$x[groups_p == "nP", 1:2], pch = "N", col = (rld.kmeans.cluster+1)[groups_p == "nP"])


