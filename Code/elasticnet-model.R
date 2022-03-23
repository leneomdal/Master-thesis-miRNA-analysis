library(glmnet)
library(caret)
library(foreach)
library(itertools)
library(boot)
source("Code//normalize-and-filter.R")

#Define data frame including ad response
log.cpm.ad = as.data.frame(t(log.cpm))
ad= metadata.df$ad
log.cpm.ad[,"AD"] = ad
#Define model matrix for glmnet model
mod.matrix = model.matrix(AD~., data = log.cpm.ad)
full.mod.matrix = model.matrix(AD~., data = log.cpm.ad)

# GLMNET tesing

#Define elastic net model using alpha = 0.5
enet.mod = glmnet(mod.matrix, y = ad, alpha = 0.5, family = "binomial", 
                  standardize = TRUE)
plot(enet.mod, xvar = "dev")

# Cross validation
set.seed(50)
cv.enet = cv.glmnet(mod.matrix, y = ad, alpha = 0.5, family = "binomial", 
                    standardize = TRUE, nfolds = 5, type.measure = "dev")
plot(cv.enet)

best.lambda = cv.enet$lambda.min
coefficients = coef(cv.enet, s = best.lambda)
final.coeffs = as.matrix(coefficients)[as.vector(coefficients) != 0,]
deviance.lambda = cv.enet$glmnet.fit$dev.ratio

dev =  (1-cv.enet$glmnet.fit$dev.ratio)*cv.enet$glmnet.fit$nulldev
dev = deviance(enet.mod)
plot(log(cv.enet$lambda), cv.enet$cvm)





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



#Function for nested CV NEEDS INPUT VECTOR OF LAMBDAS
nested.cv.alpha = function(mod.matrix, n.folds.outer, n.folds.inner, alphas, lambda.type = "lambda.1se"){
  nrow.mod.matrix = nrow(mod.matrix)
  folds.outer = sample(x = rep(1:n.folds.outer, ceiling(nrow.mod.matrix/n.folds.outer)), 
                       size = nrow.mod.matrix, replace = FALSE)
  list.foldid.inner = list()
  for(i in seq_len(n.folds.outer)){
    nrow.train.i = nrow.mod.matrix - sum(folds.outer == i)
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
                         alpha = alphas[j], foldid = as.vector(list.foldid.inner[[i]]) )
      pred = predict(cv.fit, newx = test.set, s = lambda.type, type = "response")
      
      dev = dev + find.dev(pred, as.numeric(ad[curr.test.fold]))
    }
    cv.alpha.df[j,] = c(alphas[j], dev/nrow.mod.matrix)
  }
  return(cv.alpha.df)
}
set.seed(50)
lol.df = nested.cv.alpha(mod.matrix, n.folds.outer, n.folds.inner, alphas)
View(lol.df)

# Paired bootstrap
bootstrap.elasticnet = function(log.cpm.ad, full.mod.matrix, n.boot, n.folds.outer, n.folds.inner, alphas, lambda.type = "lambda.1se"){
  boot.coeffs.df = data.frame(matrix(ncol = ncol(full.mod.matrix) - 1, nrow = 0))
  colnames(boot.coeffs.df) = colnames(mod.matrix[,-1])
  boot.enet.mods.df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(boot.enet.mods.df) = c("alpha", "lambda", "deviance")
  
  # define folds for CV of lambda after selecting alpha ( this should be equal for each bootrap sample)
  folds.lambda = sample(x = rep(1:n.folds.inner, ceiling(nrow(full.mod.matrix)/n.folds.inner)), 
                                      size = nrow(full.mod.matrix), replace = FALSE)
  for(i in 1:n.boot){
    boot.index = sample(1:nrow(log.cpm.ad), size = nrow(log.cpm.ad), replace = TRUE)
    boot.data = log.cpm.ad[boot.index, ]
    model.mat = model.matrix(AD~., data = boot.data)
    
    nested.cv.df = nested.cv.alpha(model.mat, n.folds.outer, n.folds.inner, alphas, lambda.type)
    
    
    #Find alpha with lowest mean deviance
    boot.alpha.best = nested.cv.df$alpha[nested.cv.df$deviance == min(nested.cv.df$deviance)]
    boot.cv.enet = cv.glmnet(full.mod.matrix, ad, family = "binomial", alpha = boot.alpha.best, 
                             foldid = folds.lambda)
    
    if(lambda.type == "lambda.1se"){
      boot.lamba.best = boot.cv.enet$lambda.1se
      boot.dev = boot.cv.enet$cvm[boot.cv.enet$lambda == boot.cv.enet$lambda.1se]
    }
    else if(lambda.type == "lambda.min"){
      boot.lamba.best = boot.cv.enet$lambda.min
      boot.dev = boot.cv.enet$cvm[boot.cv.enet$lambda == boot.cv.enet$lambda.min]
    }
    else{
      print("Choose either lambda.min or lambda.1se for regularization")
    }
    
    #Save data from this bootstrap sample
    boot.coeffs = as.matrix(coef(boot.cv.enet, s = boot.lamba.best))
    #Remove two first rows (intercepts)
    boot.coeffs.df[i,] = boot.coeffs[seq(3, nrow(boot.coeffs))]
    boot.enet.mods.df[i,] = c(alpha = boot.alpha.best, lambda = boot.lamba.best, deviance = boot.dev)
    print(paste("boostrap sample nr.", i))
  }
  return(list(boot.coeffs.df, boot.enet.mods.df))
}

#Define number of folds for nested CV of lambda and alpha
n.folds.inner = 5
n.folds.outer = 5

#Define vector of alphas to run cross validation
alphas =  seq(0.1, 0.9, 0.05)


set.seed(50)
list.bootstrap = bootstrap.elasticnet(log.cpm.ad, full.mod.matrix, n.boot = 5, n.folds.outer, n.folds.inner, alphas)
View(list.bootstrap[[2]])


#Save data from bootstrap
#write.csv(list.bootstrap[[2]], file = "bootstrap-models.csv", row.names = FALSE)
#write.csv(list.bootstrap[[1]], file = "bootstrap-coefficients.csv", row.names = FALSE)

