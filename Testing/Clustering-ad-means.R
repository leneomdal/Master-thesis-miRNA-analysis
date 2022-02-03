# Clustering based on atopic dermatitis

log.cpm.ad = log.cpm[,metadata.df$ad == 1]
View(log.cpm.ad)
View(metadata.df)
log.cpm.nad = log.cpm[,metadata.df$ad == 0]
means.ad = apply(log.cpm.ad, 1, mean)
means.nad = apply(log.cpm.nad, 1, mean)


log.cpm.means.ad.df = data.frame(ad = means.ad, nad = means.nad)
View(log.cpm.means.ad.df)

ad.kmeans = kmeans(log.cpm.means.ad.df, 2)
ad.clusters = ad.kmeans$cluster

expr_pca = prcomp(t(log.cpm.means.ad.df))
expr_loading = expr_pca$rotation[,1:2]
informative_loadings = rbind(
  head(expr_loading[order(abs(expr_loading[,1]), decreasing = TRUE),]),
  head(expr_loading[order(abs(expr_loading[,2]), decreasing = TRUE),])
)
biplot(x = expr_pca$x[,1:2], y= informative_loadings, scale=0)
plot(expr_pca$x[,1:2], type = "n")
points(expr_pca$x[groups_ad == "AD", 1:2], pch = "AD", col = (ad.clusters+1)[groups_ad == "AD"]) 
points(expr_pca$x[groups_ad == "nAD", 1:2], pch = "nAD", col = (ad.clusters+1)[groups_p == "nAD"])

# Scaling as suggested in book
centered.by.ad.log.cpm = log.cpm
centered.by.ad.log.cpm[,metadata.df$ad == 1] = log.cpm[,metadata.df$ad == 1] - means.ad
centered.by.ad.log.cpm[,metadata.df$ad == 0] = scaled.by.ad.log.cpm[,metadata.df$ad == 0] - means.nad
View(centered.by.ad.log.cpm)
