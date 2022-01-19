library(stats)
library(limma)
count.df = ten.days.sample.df
metadata.df = ten.days.meta.data
View(metadata.df)
ncol(count.df)
nrow(metadata.df)
View(count.df)


#Hiearchical clustering
h = hclust(dist(t(count.df)))
dist(count.df)
plot(h)
heatmap(as.matrix(count.df), scale = c( "row"))
?heatmap
count.df
heatmap(as.matrix(t(scale(t(count.df)))))


#pearson correlation

cor(count.df)
