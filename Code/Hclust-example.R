# Hierarchical clustering example
library(ggplot2)
# Simulate data
set.seed(123)
sim.norm.df = data.frame(x = rnorm(10), y = rnorm(10))
plot(sim.norm.df)

# cluster simulated data
hclust.sim = hclust(dist(sim.norm.df), method = "ward.D2")

# save plot of dendrogram
svg("C:\\Users\\Lene\\repositories\\master-thesis-latex\\figures\\hclust-example.svg")
plot(hclust.sim, main = NULL, sub = NA, xlab = "", ylab = "")
dev.off()


sim.norm.df$clust.num = 1:10
# save plot of simulated data points
svg("C:\\Users\\Lene\\repositories\\master-thesis-latex\\figures\\points-hclust.svg")
ggplot(data = sim.norm.df, aes(x = x, y = y, label = clust.num)) +
  ylim(c(-2, 2)) + xlim(-2,2) + geom_text(size = 5.5) + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
dev.off()
