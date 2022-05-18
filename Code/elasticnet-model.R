if(!require(glmnet)){
  install.packages("glmnet")
  library(glmnet)
}
if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}

print("Hola girl")

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
source("Code//enet-mod-functions.R")

#Define data frame including ad response
log.cpm.ad = as.data.frame(t(log.cpm))
ad= metadata.df$ad
log.cpm.ad[,"AD"] = ad
log.cpm.ad[,"matatopy"] = metadata.df$matatopy
log.cpm.ad[,"sex"] = metadata.df$sex
#log.cpm.ad[,"probiotic"] = metadata.df$probiotic

#Define model matrix for glmnet model
mod.matrix = model.matrix(AD~., data = log.cpm.ad)
full.mod.matrix = mod.matrix


# --------------------------REPEATED CV-------------------------------------------------
#Lambda sequence for CV
lambda.seq = exp(seq(1,-5,length=200))
#Times to repeat the CV
n.repeat = 10
#Define vector of alphas to run cross validation
alphas =  seq(0.1, 1, 0.05)
#Define number of folds in CV
n.folds = 10

# set.seed(345)
# repeated.cv.results = repeat.cv.function(mod.matrix, ad, alphas, lambda.seq, n.repeat, n.folds)
# 
# actual.rep.enet.model = glmnet(full.mod.matrix[,-1], ad, family = "binomial",
#                            alpha = repeated.cv.results$alpha, lambda = repeated.cv.results$lambda)
# 
# coeffs.actual.mod.rep = as.matrix(coef(actual.rep.enet.model))
# names.final.coeffs.rep = str_remove_all(names(as.matrix(coeffs.actual.mod.rep)[as.vector(coeffs.actual.mod.rep) != 0,]), "`")
# 

# 
# 
# # BOOTSTRAP REPEATED CV
# 
#Define number of bootstrap samples
# n.boot = 1000
# set.seed(678)
# boot.repeated.df = bootstrap.repeated.cv(log.cpm.ad, full.mod.matrix, ad, alphas,lambda.seq, n.boot, n.repeat, n.folds)
# 
# #Save data from bootstrap
# write.csv(boot.repeated.df[[2]], file = paste("Data/bootstrap-models-allcov-repeated-folds-",
#                                               n.folds,".csv",sep = ""), row.names = FALSE)
# write.csv(boot.repeated.df[[1]], file = paste("Data/bootstrap-coefficients-allcov-repeated-",
#                                               n.folds, ".csv", sep = ""), row.names = FALSE)


#----------------------NESTED CV-------------------------------------------------




#Define number of folds for nested CV of lambda and alpha
n.folds.inner = 10
n.folds.outer = 10 
lambda.type = "lambda.1se"
alphas =  seq(0.1, 1, 0.05)

#______________________

# alpha.test.nest = c()
# for(i in seq_len(50)){
#   nested.cv.alpha.df = nested.cv.alpha(full.mod.matrix, ad, n.folds.outer, n.folds.inner, alphas, lambda.type = lambda.type)
#   
#   best.alpha.nested = nested.cv.alpha.df$alpha[nested.cv.alpha.df$deviance == min(nested.cv.alpha.df$deviance)]
#   alpha.test.nest[i] = best.alpha.nested
# }
# 
# #write.csv(alpha.test.nest, file = "Data/alpha-test-nest.csv")
# test = read.csv("Data/alpha-test-nest.csv")
# View(test)
# hist(test$x)
#____________________


# RUN for actual model fit

# set.seed(2345)
# nested.cv.alpha.df = nested.cv.alpha(full.mod.matrix, ad, n.folds.outer, n.folds.inner, 
#                                      alphas, lambda.type = lambda.type)
# best.alpha.nested = nested.cv.alpha.df$alpha[nested.cv.alpha.df$deviance == min(nested.cv.alpha.df$deviance)]
# 
# final.model.cv = cv.glmnet(full.mod.matrix[,-1], ad, family = "binomial", alpha = best.alpha.nested)
# plot(final.model.cv$glmnet.fit, xvar = "lambda", label = TRUE)
# 
# 
# coeffs.final.mod = as.matrix(coef(final.model.cv, s = lambda.type))
# names.final.coeffs = str_remove_all(names(as.matrix(coeffs.final.mod)[as.vector(coeffs.final.mod) != 0,]), "`")



# BOOTSTRAP NESTED CV
set.seed(1235)
list.bootstrap = bootstrap.elasticnet(log.cpm.ad, full.mod.matrix, response = ad, n.boot = 1000,
                                      n.folds.outer, n.folds.inner, alphas,
                                      lambda.type = lambda.type)

# Define lambda type to use in name of stored bootstrap file
if(lambda.type == "lambda.1se"){
  lambda = "lambda-1se"
}
if(lambda.type == "lambda.min"){
  lambda = "lambda-min"
}

# #Save data from bootstrap
write.csv(list.bootstrap[[2]], file = paste("Data/bootstrap-models-allcov-innerF",n.folds.inner,
                                            "-outerF", n.folds.outer,"-", lambda, "-new",
                                            ".csv",sep = ""), row.names = FALSE)
write.csv(list.bootstrap[[1]], file = paste("Data/bootstrap-coefficients-allcov-innerF",
                                            n.folds.inner, "-outerF", n.folds.outer,"-",
                                            lambda,"-new", ".csv", sep = ""), row.names = FALSE)

