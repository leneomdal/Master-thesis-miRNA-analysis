if(!require(glmnet)){
  install.packages("glmnet")
  library(glmnet)
}
if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}
#library(caret)
#library(foreach)
#library(itertools)
#library(boot)
#library(tictoc)

print("Hei her er jeg! im loo cv")

project.dir = "Master-thesis-miRNA-analysis"
reg = regexpr(pattern = project.dir, getwd())
setwd(substr(getwd(), 1, reg + attr(reg, "match.length")))

path.log.cpm = "Data/log_cpm.csv"
path.metadata = "Data/metadata.csv"
if(file.exists(path.log.cpm) & file.exists(path.metadata)){
  log.cpm = read.csv(path.log.cpm)
  rownames(log.cpm) = log.cpm[,1]
  log.cpm = log.cpm[,-1]
  colnames(log.cpm) = str_remove(colnames(log.cpm), "^X")
  metadata.df = read.csv(path.metadata)
} else{
  source("Code//normalize-and-filter.R")
  write.csv(log.cpm, path.log.cpm)
  write.csv(metadata.df, path.metadata, row.names = FALSE)
}

#Define data frame including ad response
log.cpm.ad = as.data.frame(t(log.cpm))
ad= metadata.df$ad
log.cpm.ad[,"AD"] = ad
#Define model matrix for glmnet model
mod.matrix = model.matrix(AD~., data = log.cpm.ad)
full.mod.matrix = mod.matrix



# GLMNET testing
set.seed(51)
train.ind = sample(1:50, 50, replace = FALSE)
train.mat = mod.matrix[train.ind,]
test.mat = mod.matrix[c(51:60),]

#Define elastic net model using alpha = 0.5
enet.mod = glmnet(mod.matrix, y = ad, alpha = 0.5, family = "binomial", 
                  standardize = TRUE)
#plot(enet.mod, xvar = "dev")

# Cross validation
set.seed(51)
cv.enet = cv.glmnet(train.mat, y = ad[train.ind], alpha = 0.1, family = "binomial", 
                    standardize = TRUE, nfolds = 5, type.measure = "dev")

plot(cv.enet)

pred = predict(cv.enet, test.mat, s = "lambda.min",type = "response")
tru = ad[51:60]

best.lambda = cv.enet$lambda.min
coefficients = coef(cv.enet, s = best.lambda)
final.coeffs = as.matrix(coefficients)[as.vector(coefficients) != 0,]
deviance.lambda = cv.enet$glmnet.fit$dev.ratio

dev =  (1-cv.enet$glmnet.fit$dev.ratio)*cv.enet$glmnet.fit$nulldev
dev = deviance(enet.mod)
#plot(log(cv.enet$lambda), cv.enet$cvm)

# Finish testing



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

# Functions for testing speed of code
tic = function() {
  .GlobalEnv$tictoc_var = Sys.time()
}
toc = function() {
  print( Sys.time() - tictoc_var)
}


#Function for nested CV NEEDS INPUT VECTOR OF LAMBDAS
nested.cv.alpha = function(mod.matrix, response.var, n.folds.outer, n.folds.inner, alphas, lambda.type = "lambda.min"){
  nrow.mod.matrix = nrow(mod.matrix)
  folds.outer = sample(x = rep(1:n.folds.outer, ceiling(nrow.mod.matrix/n.folds.outer))[1:nrow.mod.matrix], 
                       size = nrow.mod.matrix, replace = FALSE)
  list.foldid.inner = list()
  for(i in seq_len(n.folds.outer)){
    nrow.train.i = nrow.mod.matrix - sum(folds.outer == i)
    list.foldid.inner[[i]] = sample(rep(1:n.folds.inner, 
                                        times = ceiling(nrow.train.i/n.folds.inner))[1:nrow.train.i], 
                                    nrow.train.i, replace = FALSE)
  }
  cv.alpha.df = data.frame(matrix(nrow = length(alphas), ncol = 2))
  colnames(cv.alpha.df) = c("alpha", "deviance")
  tic()
  for(j in seq_along(alphas)){                             # for each alpha run nested CV
    dev = 0
    for(i in seq_len(n.folds.outer)){                      # outer fold, nested CV 
      curr.test.fold = folds.outer == i
      curr.train.fold = folds.outer != i
      test.set = mod.matrix[curr.test.fold,]
      train.set = mod.matrix[curr.train.fold, ]
      
      cv.fit = cv.glmnet(train.set, response.var[curr.train.fold], family = "binomial",
                         alpha = alphas[j], foldid = as.vector(list.foldid.inner[[i]]),
                         standardize = TRUE)               # inner fold, CV for lambda
      
      pred = predict(cv.fit, newx = test.set, s = lambda.type, 
                     type = "response")                    # predict on left out test set of CV
      
      dev = dev + find.dev(pred, as.numeric(as.character(response.var[curr.test.fold])))
    }
    
    cv.alpha.df[j,] = c(alphas[j], dev/nrow.mod.matrix)   # store average deviance across folds
                                                          # together with its alpha value
  }
  toc()
  print("im here")
  return(cv.alpha.df)
}

#Define number of folds for nested CV of lambda and alpha
n.folds.inner = 10
n.folds.outer = nrow(full.mod.matrix)


# RUN for actual model fit

# set.seed(2345)
# nested.cv.alpha.df = nested.cv.alpha(full.mod.matrix, ad, n.folds.outer, n.folds.inner, alphas, lambda.type = "lambda.min")
# View(nested.cv.alpha.df)
# best.alpha.nested = nested.cv.alpha.df$alpha[nested.cv.alpha.df$deviance == min(nested.cv.alpha.df$deviance)]
# 
# final.model.cv = cv.glmnet(full.mod.matrix, ad, family = "binomial", alpha = best.alpha.nested)
# plot(final.model.cv$glmnet.fit, xvar = "lambda", label = TRUE)
# final.model.cv$lambda.1se
# 
coeffs.final.mod = as.matrix(coef(final.model.cv, s = "lambda.1se"))
names.final.coeffs = str_remove_all(names(as.matrix(coeffs.final.mod)[as.vector(coeffs.final.mod) != 0,]), "`")




# Paired bootstrap
bootstrap.elasticnet = function(log.cpm.ad, full.mod.matrix, n.boot, n.folds.outer, 
                                n.folds.inner, alphas, lambda.type = "lambda.min"){
  boot.coeffs.df = data.frame(matrix(ncol = ncol(full.mod.matrix) - 1, nrow = 0))
  colnames(boot.coeffs.df) = colnames(mod.matrix[,-1])
  boot.enet.mods.df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(boot.enet.mods.df) = c("alpha", "lambda", "deviance")
  
  # define folds for CV of lambda after selecting alpha ( this should be equal for each bootstrap sample)
  folds.lambda = sample(x = rep(1:n.folds.inner, ceiling(nrow(full.mod.matrix)/n.folds.inner)), 
                                      size = nrow(full.mod.matrix), replace = FALSE)
  for(i in 1:n.boot){
    boot.index = sample(1:nrow(log.cpm.ad), size = nrow(log.cpm.ad), replace = TRUE)
    boot.data = log.cpm.ad[boot.index, ]
    boot.mod.matrix = model.matrix(AD~., data = boot.data)
    boot.ad = ad[boot.index]
    
    nested.cv.df = nested.cv.alpha(boot.mod.matrix, boot.ad, n.folds.outer, n.folds.inner, 
                                   alphas, lambda.type = lambda.type)
    
    #Find alpha with lowest mean deviance
    boot.alpha.best = nested.cv.df$alpha[nested.cv.df$deviance == min(nested.cv.df$deviance)]
    #If multiple alphas show equal deviance, choose the highest alpha
    if(lenght(boot.alpha.best) > 1){
      boot.alpha.best = boot.alpha.best[length(boot.alpha.best)]
    }
    
    
    # SKAL MAN BRUKE FULL.MOD.MAT HER ELLER HELE BOOTSTRAPPED DATASET??
    boot.cv.enet = cv.glmnet(boot.mod.matrix, boot.ad, family = "binomial", alpha = boot.alpha.best, 
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
    #Remove two first rows (intercepts) DO NOT DO THIS
    #boot.coeffs.df[i,] = boot.coeffs[seq(3, nrow(boot.coeffs))]
    boot.enet.mods.df[i,] = c(alpha = boot.alpha.best, lambda = boot.lamba.best, 
                              deviance = boot.dev)
    print(paste("boostrap sample nr.", i))
  }
  return(list(boot.coeffs.df, boot.enet.mods.df))
}



#Define vector of alphas to run cross validation
alphas =  seq(0.1, 1, 0.05)
lambda.t = "lambda.1se"

set.seed(1235)
list.bootstrap = bootstrap.elasticnet(log.cpm.ad, full.mod.matrix, n.boot = 1000, 
                                      n.folds.outer, n.folds.inner, alphas, 
                                      lambda.type = lambda.t)




#Save data from bootstrap
write.csv(list.bootstrap[[2]], file = paste("Data/bootstrap-models-innerF",n.folds.inner, 
                                            "-outerF", n.folds.outer,"-", lambda.t, 
                                            ".csv",sep = ""), row.names = FALSE)
write.csv(list.bootstrap[[1]], file = paste("Data/bootstrap-coefficients-innerF",
                                            n.folds.inner, "-outerF", n.folds.outer,"-", 
                                            lambda.t ,".csv", sep = ""), row.names = FALSE)

