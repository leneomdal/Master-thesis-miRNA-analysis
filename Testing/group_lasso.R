#group lasso
library(gglasso)
library(tictoc)
source("Code//new-heatmaps.R")
plot(clust.m.cor)
as.dendrogram(clust.m.cor) %>% plot()
clustering.groups.cor = cutree(tree = as.dendrogram(clust.m.cor), k = 8)
clustering.groups.cor

#Define data frame including ad response
log.cpm.ad.gg = as.data.frame(t(log.cpm))
ad= metadata.df$ad
log.cpm.ad.gg[,"AD"] = ad
log.cpm.ad.gg[,"matatopy"] = as.factor(ifelse(metadata.df$matatopy == 0, 0, 1))
log.cpm.ad.gg[,"noMatatopy"] = as.factor(ifelse(metadata.df$matatopy == 1, 0, 1))
log.cpm.ad.gg[,"boy"] = as.factor(ifelse(metadata.df$sex == 0, 0, 1))
log.cpm.ad.gg[,"girl"] = as.factor(ifelse(metadata.df$sex == 1, 0, 1))

#log.cpm.ad.gg[,"probiotic"] = metadata.df$probiotic

#Define model matrix for glmnet model
mod.matrix.gg = model.matrix(AD~., data = log.cpm.ad.gg)

ad.gglasso = ifelse(ad == 0, -1, 1)
#ad.gglasso = as.factor(ad)


group.lasso.fit = gglasso(mod.matrix.gg[,-1], ad.gglasso, 
                          group = c(clustering.groups.cor, 9,9, 10, 10), 
                          loss = "logit", intercept = TRUE)
plot(group.lasso.fit)
tic()
cv.group.lasso.fit = cv.gglasso(mod.matrix.gg[,-1], ad.gglasso, 
                                group = c(clustering.groups.cor, 9,9, 10, 10), 
                                loss = "logit", intercept = TRUE)
toc()
plot(cv.group.lasso.fit)

cv.group.lasso.fit$lambda.min
cv.group.lasso.fit$lambda.1se

#Find coefficients

coeffs.group.lasso.min =as.matrix(coef(cv.group.lasso.fit$gglasso.fit, 
                                       s = cv.group.lasso.fit$lambda.min))
coeffs.group.lasso.1se =as.matrix(coef(cv.group.lasso.fit$gglasso.fit, 
                                       s = cv.group.lasso.fit$lambda.1se))

#save models
write.csv(coeffs.group.lasso.min, file = "Data/group-lasso-min-cor-groups-8-min.csv", 
          row.names = FALSE)
write.csv(coeffs.group.lasso.1se, file = "Data/group-lasso-min-cor-groups-8-1se.csv", 
          row.names = FALSE)


rownames(coeffs.group.lasso.1se)[coeffs.group.lasso.min != 0]


