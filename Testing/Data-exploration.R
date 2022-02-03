plot(rowMeans(count.df))
plot(rowMeans(cpm))
plot(rowMeans(log.cpm))
centered.log.cpm = t(scale(t(log.cpm), center = TRUE, scale = FALSE))
scaled.log.cpm = t(scale(t(log.cpm), center = TRUE, scale = TRUE))

plot(rowMeans(centered.log.cpm))
plot(rowMeans(scaled.log.cpm))

#Heatmap euclidian distance with scaled log cpm values
create.heatmap(scaled.log.cpm, type = "miRNA",)
