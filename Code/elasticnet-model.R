library(glmnet)
library(caret)
library(foreach)
library(itertools)
library(boot)

# GLMNET
log.cpm.ad = as.data.frame(t(log.cpm))
ad= metadata.df$ad
log.cpm.ad[,"AD"] = ad
mod.matrix = model.matrix(AD~., data = log.cpm.ad)

#Define elsatic net model using alpha = 0.5
enet.mod = glmnet(mod.matrix, y = ad, alpha = 0.5, family = "binomial", 
                  standardize = TRUE)
plot(enet.mod, xvar = "dev")
enet.mod$lambda


# Cross validation
set.seed(50)
cv.enet = cv.glmnet(mod.matrix, y = ad, alpha = 0.5, family = "binomial", 
                    standardize = TRUE, nfolds = 5, type.measure = "dev")
plot(cv.enet)
cv.enet$cvm


best.lambda = cv.enet$lambda.min
coefficients = coef(cv.enet, s = best.lambda)
final.coeffs = as.matrix(coefficients)[as.vector(coefficients) != 0,]
deviance.lambda = cv.enet$glmnet.fit$dev.ratio

dev =  (1-cv.enet$glmnet.fit$dev.ratio)*cv.enet$glmnet.fit$nulldev
dev = deviance(enet.mod)
plot(log(cv.enet$lambda), cv.enet$cvm)

#Define folds
set.seed(50)
#Define vector of alphas to run cross validation
alphas =  seq(0.1, 0.9, 0.05)
#Define vector of lambdas
lambdas = seq(0,1)




# cross validate lambdas  NEEED TO DEFINE A VECTOR OF LAMBDAS HERE!!!! DEFINE FOLDID, SAME FOLDS FOR EACH ALPHA
# cv.df = data.frame(alpha = rep(0, length(alphas)), mcve = rep(0,length(alphas)), lambda = rep(0,length(alphas)))
# View(cv.df)
# for(i in 1:length(alphas)){
#   cv = cv.glmnet(mod.matrix, ad, type.measure = "deviance", family = "binomial", 
#                  standardize = TRUE, alpha = alphas[i], foldid = folds)
#   cv.df[i,] = c(alphas[i], cv$cvm[cv$lambda ==cv$lambda.min], cv$lambda.min)
# }
# 
# alpha.min = cv.df$alpha[cv.df$mcve == min(cv.df$mcve)]
# lambda.min = cv.df$lambda[cv.df$mcve == min(cv.df$mcve)] 
# 
# enet.fit = glmnet(mod.matrix, ad, family = "binomial", standardize = TRUE, 
#                   lambda = lambda.min, alpha = alpha.min)
# enet.fit$lambda
# coef(enet.fit)


# Function for calculating deviance
find.dev = function(pred, truth){
  dev = -2*(truth %*% log(pred) + (1 - truth) %*% log(1 - pred))
  return(dev)
}


#Define number of folds for nested CV of lambda and alpha
n.folds = 5
#Function for nested cv NEEDS INPUT VECTOR OF LAMBDAS
nested.cv.alpha = function(n.folds, alphas, lambda.type = "lambda.1se"){
  folds.alpha = sample(rep(1:n.folds, 12), nrow(mod.matrix), replace = FALSE)
  folds.lambda = sample(rep(1:n.folds, 12), nrow(mod.matrix) - (nrow(mod.matrix)/n.folds), 
                        replace = FALSE)
  cv.alpha.df = data.frame(alpha = rep(NA,length(alphas)), deviance = rep(NA, length(alphas)))
  
  for(j in seq_along(alphas)){
    dev = 0
    for(i in 1:n.folds){
      curr.test.fold = folds.alpha == i
      curr.train.fold = folds.alpha != i
      test.set = mod.matrix[curr.test.fold,]
      train.set = mod.matrix[curr.train.fold, ]
      cv.fit = cv.glmnet(train.set, ad[curr.train.fold], family = "binomial", 
                         standadize = TRUE, alpha = a, foldid = folds.lambda )
      if(lambda.type == "lambda.1se"){
        lambda = cv.fit$lambda.1se
      }
      else if(lambda.type == "lambda.min"){
        lambda = cv.fit$lambda.min
      }
      else{
        print("Choose either lambda.min or lambda.1se for regularization")
      }
      best.fit = glmnet(train.set, ad[curr.train.fold], family = "binomial", alpha = alphas[j],
                        lambda = lambda, foldid = folds.lambda)
      pred = predict(best.fit, test.set, s = lambda, type = "response")
      dev = dev + find.dev(pred, ad[curr.test.fold])
    }
    cv.alpha.df[j,] = c(alphas[j], dev/nrow(mod.matrix))
  }
  return(cv.alpha.df)
}

set.seed(50)
nested.cv.alpha(n.folds, alphas, lambda.type = "lambda.1se")

#Find alpha min from nested cv
alpha.min = cv.alpha.df$alpha[cv.alpha.df$deviance == min(cv.alpha.df$deviance)]
# Cross validate for lambdas using alpha min on entire data set
cv.enet = cv.glmnet(mod.matrix, ad, family = "binomial", alpha = alpha.min)
plot(cv.enet)

lambda.se = cv.enet$lambda.1se

enet.mod = glmnet(mod.matrix, ad, family = "binomial", alpha = alpha.min, 
                  lambda = lambda.se)
coeffs = coef(enet.mod)
included.coeffs = as.matrix(coeffs)[as.vector(coeffs) != 0,]


# Paired bootstrap
?boot

bootstrap.elasticnet = function(n.boot, alpha, lambda){
  included.coeffs =c()
  for(i in 1:n.boot){
    boot.index = sample(1:ncol(log.cpm), size = ncol(log.cpm), replace = TRUE)
    boot.df = log.cpm.ad[boot.index, ]
    
    model.mat = model.matrix(AD~., data = boot.df)
    elasticnet.mod.boot = glmnet(model.mat, ad[boot.index], family = "binomial", 
                                 alpha = alpha, lambda = lambda)
    coeffs = coef(elasticnet.mod.boot)
    included.coeffs = c(as.matrix(coeffs)[as.vector(coeffs) != 0,][-1], included.coeffs)
    curr.coeffs = as.matrix(coeffs)[as.vector(coeffs) != 0,][-1]
    for(j in seq_along(curr.coeffs)){
      
    }
  }
  return(included.coeffs)
}
coeffs.boot = bootstrap.elasticnet(100, alpha.min, lambda.se)



hist(coeffs.boot[names(coeffs.boot) == "`miR-20a-3p`"])
coeffs.boot

matrix(data= c(1,2,3,4,5,6), nrow=2, ncol = 3)[]
