source("Code//normalize-and-filter.R")
library(svglite)
#source("Testing//hierarchical-clustering-functions.R")
plot(rowMeans(count.df))
plot(rowMeans(cpm))
plot(rowMeans(log.cpm))
centered.log.cpm = t(scale(t(log.cpm), center = TRUE, scale = FALSE))
scaled.log.cpm = t(scale(t(log.cpm), center = TRUE, scale = TRUE))

plot(rowMeans(centered.log.cpm))
plot(rowMeans(scaled.log.cpm))



# Histogram of cpm values for miRNA  7 and 84
df.mirna = data.frame(mirna7 = cpm[rownames(cpm) == "let-7d-3p"], mirna80 = cpm[rownames(cpm) == "miR-135a-5p"])
rownames(cpm)[84]

#save plots
svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\hist-cpm-mirna7.svg")
ggplot(df.mirna, aes(x= mirna7)) + geom_histogram(binwidth=40, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of let-7d-3p cpm values") + xlab("cpm values") + ylab("Count") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#Close the graphics device
dev.off() 

svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\hist-cpm-mirna84.svg")
ggplot(df.mirna, aes(x= mirna80)) + geom_histogram(binwidth=25, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of miR-135a-5p cpm values") + xlab("cpm values") + ylab("Count") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


# Histogram of log cpm values for miRNA  7 and 84
df.log.mirna = data.frame(mirna7 = log.cpm[rownames(log.cpm) == "let-7d-3p"], mirna80 = log.cpm[rownames(log.cpm) == "miR-135a-5p"])

svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\hist-log-cpm-mirna7.svg")
ggplot(df.log.mirna, aes(x= mirna7)) + geom_histogram(binwidth=0.08, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of let-7d-3p log cpm values") + xlab("log cpm values") + ylab("Count") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\hist-log-cpm-mirna84.svg")
ggplot(df.log.mirna, aes(x= mirna80)) + geom_histogram(binwidth=0.14, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of miR-135a-5p log cpm values") + xlab("log cpm values") + ylab("Count") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



# HISTOGRAM OF MEAN VALUES OF MIRNA

mean.mirnas.df = data.frame(means = rowMeans(log.cpm))
svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\hist-mean-mirna.svg")
ggplot(mean.mirnas.df, aes(x = means)) + geom_histogram(binwidth=0.3, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of sample means for each miRNA") + xlab("log cpm mean") + ylab("Count") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# HISTOGRAM OF VARIANCES OF MIRNA
var.mirnas.df = data.frame(var = apply(log.cpm, 1, var))

#mirnas with highest variance
rownames(var.mirnas.df)[var.mirnas.df == max(var.mirnas.df$var)]
rownames(var.mirnas.df)[order(var.mirnas.df$var, decreasing = TRUE)[1:5]]


mean.var.mirnas = mean(var.mirnas.df$var)

estimated.degrees.of.freedom = mean.var.mirnas/2


# PLOT EMIRICAL VARIANCE
svg("C:\\Users\\Lene\\Documents\\Skole\\Master\\Master-thesis-miRNA-analysis\\Figures\\histogram-var-mirna.svg")
ggplot(var.mirnas.df, aes(x = var)) + geom_histogram(binwidth=0.1, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  ggtitle("Histogram of sample variance for each miRNA") + xlab("Variance") + ylab("Count") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# PLOT POSTERIOR S
hist(fit.vanilla$s2.post)
ggplot() + geom_histogram(data = data.frame(post = fit.vanilla$s2.post), aes(x = post), binwidth=0.1, fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
                            ggtitle("Histogram of var miRNA values") + xlab("var") + ylab("Count") + geom_vline(xintercept = s_0) 



# GAMMA DISTRIBUTION
x_lower_g <- 0
x_upper_g <- 5

# Gamma Distribution Plot With Rate = 2 and Scale = 0.5

ggplot(data.frame(x = c(x_lower_g , x_upper_g)), aes(x = x)) + 
  xlim(c(x_lower_g , x_upper_g)) + 
  stat_function(fun = dgamma, args = list(rate = 1/2, shape = estimated.degrees.of.freedom/2), geom = "area", 
                fill = "orange", alpha = 0.25) + 
  stat_function(fun = dgamma, args = list(rate = 1/2, shape = estimated.degrees.of.freedom/2)) + 
  labs(x = "\n x", y = "f(x) \n", 
       title = "Gamma Distribution With Rate & Shape = 1/2 \n")

curve(1/(d_0*s_0^2)*dchisq(x, df = d_0), from = 0, to = 40)

# plot of variance against mean values

plotSA(eB.voom.fit, main="Mean variance points and prior variance")
?plotSA


# correlation plot

library(corrplot)
library(RColorBrewer)
cor.mirna = cor(scale(t(log.cpm[1:10,])), method = "spearman")

corrplot(cor.mirna)

# Colors
corrplot(cor.mirna, method = "color",
         title = "method = 'color'",
         tl.pos = "n", mar = c(2, 1, 3, 1)) 


corrplot(cor.mirna,
         method = "circle",       
         order = "hclust",         # Ordering method of the matrix
         hclust.method = "ward.D", # If order = "hclust", is the cluster method to be used
         addrect = 2,              # If order = "hclust", number of cluster rectangles
         rect.col = 3,             # Color of the rectangles
         rect.lwd = 3)             # Line width of the rectangles




# BOXPLOT OF MOST ABUNDANT miRNAs
sum.reads = apply(dge$counts, 1, sum)
sorted.sum.reads = sum.reads[order(sum.reads, decreasing = TRUE)]
n.mirna.to.include = 20
most.abundant.mirna.names = names(sorted.sum.reads[seq(1, n.mirna.to.include)])
most.abundant.mirnas = dge$counts[most.abundant.mirna.names,]
abundant.mirna.long = pivot_longer(data = as.data.frame(t(most.abundant.mirnas)), cols = rownames(most.abundant.mirnas), names_to = "miRNA", values_to = "cpmValue")
#abundant.mirna.long$miRNA = str_remove(abundant.mirna.long$miRNA, "hsa-")
abundant.mirna.long$miRNA = factor(abundant.mirna.long$miRNA, levels = most.abundant.mirna.names)

ggplot(data = abundant.mirna.long, aes(x = miRNA, y = cpmValue )) + geom_boxplot() + ylab("cpm value") + ggtitle("20 most abundant miRNAs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) 
ggsave("C:\\Users\\Lene\\repositories\\master-thesis-latex\\figures\\most-abundant.svg", height = 4.5, width = 5)


# SEPARATED BOXPLOT

abundant.mirna.long$probiotic.bool = rep(ten.days.meta.data$probiotic, each = n.mirna.to.include)

abundant.mirna.long$Probiotic = ifelse(abundant.mirna.long$probiotic.bool == 1, "Yes", "No" )




ggplot(data = abundant.mirna.long, aes(x = miRNA, y = cpmValue, fill = Probiotic)) + geom_boxplot() + 
  ylab("cpm value") + ggtitle("20 most abundant miRNAs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.major.x = element_blank())
ggsave("C:\\Users\\Lene\\repositories\\master-thesis-latex\\figures\\most-abundant-sep.svg", height = 4.5, width = 6)



# LIBRARY SIZES 

lib.sizes = apply(dge$counts, 2, sum)

svg("C:\\Users\\Lene\\repositories\\master-thesis-latex\\figures\\lib-sizes.svg")
ggplot(data = data.frame(lib = lib.sizes, mirna = names(lib.sizes)), aes(x = mirna, y = lib)) + 
  geom_bar(stat = "identity",  fill="#69b3a2", alpha=0.9) + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.border = element_blank()) +
  xlab("Samples") + ylab("Library size") + ggtitle("Bar plot of library sizes") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
  

(metadata.df$sens2yr == 1 & metadata.df$probiotic == 0)

