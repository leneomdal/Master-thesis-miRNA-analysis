library(glmnet)
library(elasticnet)



# GLMNET
log.cpm.df = as.data.frame(t(log.cpm))
y = metadata.df$ad
log.cpm.df[,"AD"] = y


cpm.model.matrix = model.matrix(AD~., data = log.cpm.df)

View(cpm.model.matrix)


elasticnet.mod = glmnet(cpm.model.matrix, y, alpha = 0.8, family = "binomial", standardize = TRUE)
plot(elasticnet.mod)


# CV GLMNET
set.seed(1)
cv.out = cv.glmnet(cpm.model.matrix, y , alpha=1, family = "binomial", standardize = TRUE)
plot(cv.out)

best.lambda.lasso = cv.out$lambda.min
cv.out$glmnet.fit$lambda
best.lambda.pos = which(cv.out$glmnet.fit$lambda == best.lambda.lasso)

which(cv.out$glmnet.fit$beta[,15] != 0)

sum(coef(cv.out, s = best.lambda.lasso) != 0)
