library(glmnet)
library(caret)
library(foreach)
library(itertools)

# GLMNET
log.cpm.df = as.data.frame(t(log.cpm))
ad= metadata.df$ad
log.cpm.df[,"AD"] = ad
mod.matrix = model.matrix(AD~., data = log.cpm.df)

#Define elsatic net model using alpha = 0.5
enet.mod = glmnet(mod.matrix, y = ad, alpha = 0.5, family = "binomial", 
                  standardize = TRUE)
plot(enet.mod, xvar = "dev")
enet.mod$dev.ratio


# Cross validation
set.seed(50)
cv.enet = cv.glmnet(mod.matrix, y = ad, alpha = 0.5, family = "binomial", 
                    standardize = TRUE, nfolds = 5, type.measure = "dev")
plot(cv.enet)
?plot.cv.glmnet


best.lambda = cv.enet$lambda.min
coefficients = coef(cv.enet, s = best.lambda)
final.coeffs = as.matrix(coefficients)[as.vector(coefficients) != 0,]
deviance.lambda = cv.enet$glmnet.fit$dev.ratio

dev =  (1-cv.enet$glmnet.fit$dev.ratio)*cv.enet$glmnet.fit$nulldev
dev = deviance(enet.mod)
plot(log(cv.enet$lambda), dev)

#Define folds
folds = sample(1:5,size=length(ad),replace=TRUE)
#Define vector of alphas to run cross validation
alphas = seq(0, 1, 10)
#Define vector og lambdas
lambdas = 


train_control <- trainControl(method = "repeatedcv",
                                number = 5,
                                repeats = 5,
                                search = "random",
                                verboseIter = TRUE)
?trainControl

train_control
# Train the model
elastic_net_model <- train(mpg ~ .,
                           data = cbind(y, X),
                           method = "glmnet",
                           preProcess = c("center", "scale"),
                           tuneLength = 25,
                           trControl = train_control)


# cross validate alpha
alphas =  seq(0.1, 0.9, 0.05)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(mod.matrix, ad, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
}

cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(mod.matrix, ad, family = "binomial", lambda = cv3$lambda.min, alpha = cv3$alpha)
coef(md3)




# cross validate lambdas
cv.df = data.frame(alpha = rep(0, length(alphas)), mcve = rep(0,length(alphas)), lambda = rep(0,length(alphas)))
View(cv.df)
for(i in 1:length(alphas)){
  cv = cv.glmnet(mod.matrix, ad, type.measure = "deviance", family = "binomial", 
                 standardize = TRUE, alpha = alphas[i])
  cv.df[i,] = c(alphas[i], cv$cvm[cv$lambda ==cv$lambda.min], cv$lambda.min)
}

alpha.min = cv.df$alpha[cv.df$mcve == min(cv.df$mcve)]
lambda.min = cv.df$lambda[cv.df$mcve == min(cv.df$mcve)] 




