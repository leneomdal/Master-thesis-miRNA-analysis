source("Code//normalize-and-filter.R")
#source("Testing//hierarchical-clustering-functions.R")
plot(rowMeans(count.df))
plot(rowMeans(cpm))
plot(rowMeans(log.cpm))
centered.log.cpm = t(scale(t(log.cpm), center = TRUE, scale = FALSE))
scaled.log.cpm = t(scale(t(log.cpm), center = TRUE, scale = TRUE))

plot(rowMeans(centered.log.cpm))
plot(rowMeans(scaled.log.cpm))






# Histogram of cpm values for miRNA  7 and 84
df.mirna = data.frame(mirna7 = cpm[7,], mirna80 = cpm[84,])
rownames(cpm)[84]

#save plots
svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\histogram-cpm-mirna7.svg")
ggplot(df.mirna, aes(x= mirna7)) + geom_histogram(binwidth=40, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of let-7d-3p cpm values") + xlab("cpm values") + ylab("Count")
#Close the graphics device
dev.off() 

svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\histogram-cpm-mirna84.svg")
ggplot(df.mirna, aes(x= mirna80)) + geom_histogram(binwidth=40, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of miR-135a-5p cpm values") + xlab("cpm values") + ylab("Count")
dev.off()


# Histogram of log cpm values for miRNA  7 and 84
df.log.mirna = data.frame(mirna7 = log.cpm[7,], mirna80 = log.cpm[84,])

svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\histogram-log-cpm-mirna7.svg")
ggplot(df.log.mirna, aes(x= mirna7)) + geom_histogram(binwidth=0.1, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of let-7d-3p log cpm values") + xlab("log cpm values") + ylab("Count")
dev.off()

svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\histogram-log-cpm-mirna84.svg")
ggplot(df.log.mirna, aes(x= mirna80)) + geom_histogram(binwidth=0.2, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of miR-135a-5p log cpm values") + xlab("log cpm values") + ylab("Count")
dev.off()

