source("Code//new-heatmaps.R")
source("Code//heatmap-top-mirna.R")
library(xtable)

# Exploring the clusters

#plot dendrograms
as.dendrogram(clust.s.euclid) %>% plot()
as.dendrogram(clust.s.cor) %>% plot()
as.dendrogram(clust.s.abs.cor) %>% plot()
as.dendrogram(clust.s.top) %>% plot()
plot(clust.s.abs.cor)


# Function for extracting intormation about the samples in each cluster
get.table.summary.clust = function(hclust.obj, k){
  clusters = cutree(tree = as.dendrogram(hclust.obj), k = k)
  table.df = data.frame(probiotic = rep(0, k))
  for(i in 1:k){
    table.df$probiotic[i] = sum(metadata.df$probiotic[clusters==i] == 1)
    table.df$placebo[i] = sum(metadata.df$probiotic[clusters==i] == 0)
    table.df$m.atopy[i] = sum(metadata.df$matatopy[clusters==i] == 1)
    table.df$no.m.atopy[i] = sum(metadata.df$matatopy[clusters==i] == 0)
    table.df$ad[i] = sum(metadata.df$ad[clusters==i] == 1)
    table.df$no.ad[i] = sum(metadata.df$ad[clusters==i] == 0)
    
  }
  return(table.df)
}
table.euclid = get.table.summary.clust(clust.s.euclid, k = 4)
table.cor = get.table.summary.clust(clust.s.cor, k = 5)
table.abs.cor = get.table.summary.clust(clust.s.abs.cor, k = 5)
table.top = get.table.summary.clust(clust.s.top, k = 3)

# Generate table for latex
#xtable(table.euclid)
#xtable(table.cor)
#xtable(table.abs.cor)
#xtable(table.top)



# compare clusters of euclid and cor distance
mirna.clusters.cor = cutree(tree = as.dendrogram(clust.m.cor), k = 2)
mirna.cluster.euclid = cutree(tree = as.dendrogram(clust.m.euclid), k = 2)
cluster.2.mirna.cor = which(mirna.clusters.cor == 2) 
cluster.2.mirna.euclid = which(mirna.cluster.euclid == 2)
difference = 0
for(elm in cluster.2.mirna.euclid){
  if(elm %in% cluster.2.mirna.cor){}
  else{
    difference = difference +1
    print(elm)
    
  }
}
difference
length(cluster.2.mirna.cor)
length(cluster.2.mirna.euclid)
