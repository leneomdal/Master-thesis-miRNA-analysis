library(glmnet)
library(caret)
library(foreach)
library(itertools)
library(boot)
source("Code//normalize-and-filter.R")
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
n.folds.inner = 5
n.folds.outer = 5
#Function for nested CV NEEDS INPUT VECTOR OF LAMBDAS
nested.cv.alpha = function(mod.matrix, n.folds.outer, n.folds.inner, alphas, lambda.type = "lambda.1se"){
  folds.outer = sample(x = rep(1:n.folds.outer, ceiling(nrow(mod.matrix)/n.folds.outer)), 
                       size = nrow(mod.matrix), replace = FALSE)
  list.foldid.inner = list()
  for(i in seq_len(n.folds.outer)){
    nrow.train.i = nrow(mod.matrix) - sum(folds.outer == i)
    list.foldid.inner[[i]] = sample(rep(1:n.folds.inner, times = ceiling(nrow.train.i/n.folds.inner)), 
                                    nrow.train.i, replace = FALSE)
  }
  cv.alpha.df = data.frame(alpha = rep(NA,length(alphas)), deviance = rep(NA, length(alphas)))
  
  for(j in seq_along(alphas)){
    dev = 0
    for(i in seq_len(n.folds.outer)){
      curr.test.fold = folds.outer == i
      curr.train.fold = folds.outer != i
      test.set = mod.matrix[curr.test.fold,]
      train.set = mod.matrix[curr.train.fold, ]
      cv.fit = cv.glmnet(train.set, ad[curr.train.fold], family = "binomial",
                         alpha = a, foldid = as.vector(list.foldid.inner[[i]]) )
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
                        lambda = lambda)
      pred = predict(best.fit, test.set, s = lambda, type = "response")
      
      
      dev = dev + find.dev(pred, as.numeric(ad[curr.test.fold]))
    }
    cv.alpha.df[j,] = c(alphas[j], dev/nrow(mod.matrix))
  }
  return(cv.alpha.df)
}

set.seed(50)
nested.cv.alpha(mod.matrix, n.folds.outer, n.folds.inner, alphas, lambda.type = "lambda.1se")

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
bootstrap.elasticnet = function(log.cpm.ad, n.boot,){
  coeffs.boot.df = data.frame( row.names = rownames(log.cpm))
  
  for(i in 1:n.boot){
    boot.index = sample(1:nrow(log.cpm.ad), size = nrow(log.cpm.ad), replace = TRUE)
    boot.data = log.cpm.ad[boot.index, ]
    
    model.mat = model.matrix(AD~., data = boot.df)
    
    nested.cv.alpha()
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


coeffs.boot = bootstrap.elasticnet(1000, alpha.min, lambda.se)
table.boot = table(names(coeffs.boot))
barplot(table.boot)
hist(coeffs.boot[names(coeffs.boot) == "`miR-20a-3p`"])
