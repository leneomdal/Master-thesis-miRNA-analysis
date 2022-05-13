# group lasso
library(gglasso)

plot(clust.m.cor)
as.dendrogram(clust.m.cor) %>% plot()
clustering.groups.cor = cutree(tree = as.dendrogram(clust.m.cor), k = 8)
clustering.groups.cor

ad.gglasso = ifelse(ad == 0, -1, 1)

group.lasso.fit = gglasso(mod.matrix[,-1], ad.gglasso, group = clustering.groups.cor, loss = "logit", intercept = TRUE)
plot(group.lasso.fit)
tic()
cv.group.lasso.fit = cv.gglasso(mod.matrix[,-1], ad.gglasso, group = clustering.groups.cor, loss = "logit", intercept = TRUE)
toc()
plot(cv.group.lasso.fit)

cv.group.lasso.fit$lambda.min

all.coeffs.group.lasso =as.matrix(coef(cv.group.lasso.fit$gglasso.fit, s = cv.group.lasso.fit$lambda.min))

rownames(all.coeffs.group.lasso)[all.coeffs.group.lasso != 0]


library(polyester)

#simulate dataset
