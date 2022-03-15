rm(list=ls()) 
source("Code//load_count_data.R")
source("Code//normalize-and-filter.R")
library(biglasso)
library(glmnet)

# BIGLASSO
X.cpm = as.big.matrix(t(log.cpm))    # big matrix without intercept
X.rld = as.big.matrix(t(assay(rld)))
Y.ad = groups_ad

lasso.fit = biglasso( X.cpm, Y.ad, family = "binomial", penalty = "lasso", alg.logistic = "Newton")
plot(lasso.fit)
#plot(lasso.fit,  main = "lasso",  xlim = c(-1.5,-4.7), ylim = c(-0.5, 1.5))



# GLMNET lasso
cpm.filtered.with.AD = as.data.frame(cpm.filtered)
cpm.filtered.with.AD["AD",] = groups_ad
full.mat = as.data.frame(t(cpm.filtered.with.AD))
full.mat[,colnames(full.mat) != "AD"] = sapply(full.mat[,colnames(full.mat) != "AD"], as.numeric)
x = model.matrix(AD ~., data = full.mat)


lasso.mod = glmnet(x, full.mat$AD, alpha = 1, family = "binomial", standardize = TRUE, 
                   intercept = TRUE)
plot(lasso.mod, xvar = "lambda", xlim = c(-1.5,-4.7), ylim = c(-0.5, 1.5))


# CV GLMNET

set.seed(1)
cv.out=cv.glmnet(x, metadata.df$ad, alpha=1, family = binomial)
plot(cv.out)

best.lambda.lasso = cv.out$lambda.min



# CV biglasso
cv.fit.lcpm.lasso = cv.biglasso(X.cpm, as.numeric(as.character(metadata.df$ad)), family = "binomial", 
                     seed = 1234, nfolds = 10, ncores = 8)
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1) ,mgp = c(2.5, 0.5, 0))
plot(cv.fit.lcpm.lasso, type = "all")

cvfit <- cv.biglasso(X.cpm, as.numeric(as.character(metadata.df$ad)), family = "binomial", 
                     penalty = "enet", alpha = 0.5, seed = 1234, nfolds = 10, ncores = 4)
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1), mgp = c(2.5, 0.5, 0))
plot(cvfit, type = "all")




cv.fit.rld = cv.biglasso(X.rld, as.numeric(as.character(metadata.df$ad)), family = "binomial", 
                         penalty = "lasso", seed = 1234, nfolds = 10, ncores = 4)
par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1) ,mgp = c(2.5, 0.5, 0))
plot(cv.fit.rld, type = "all")


# Extracting linear combination from lasso log cpm
t(log.cpm)
cv.fit.lcpm.lasso$min
best.lc.lasso.lcpm = cv.fit.lcpm.lasso$fit$beta[, cv.fit.lcpm.lasso$min]
best.lc.lasso.lcpm
View(log.cpm)
which(best.lc.lasso.lcpm != 0)

sum(names(best.lc.lasso.lcpm[-1]) != row.names(log.cpm))

log.cpm.df = as.data.frame(t(log.cpm))
log.cpm.df$best.lc = best.lc.lasso.lcpm[1] + as.matrix(log.cpm.df) %*% best.lc.lasso.lcpm[-1]
log.cpm.df$best.lc
metadata.and.best.lc.df = metadata.df
metadata.and.best.lc.df$best.lc = log.cpm.df$best.lc
View(metadata.and.best.lc.df)
summary(lm(best.lc ~ matatopy + probiotic + sex + sib , metadata.and.best.lc.df))
