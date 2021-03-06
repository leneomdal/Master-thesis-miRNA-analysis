if(!require(glmnet)){
  install.packages("glmnet")
  library(glmnet)
}
if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}



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
#log.cpm.ad[,"matatopy"] = as.factor(ifelse(metadata.df$matatopy == 0, -1, 1))
#log.cpm.ad[,"sex"] = as.factor(ifelse(metadata.df$sex == 0, -1, 1))
#log.cpm.ad[,"probiotic"] = metadata.df$probiotic

#Define model matrix for glmnet model
mod.matrix = model.matrix(AD~., data = log.cpm.ad)
full.mod.matrix = mod.matrix

log.cpm.extended = log.cpm.ad
log.cpm.extended[,"matatopy"] = as.factor(ifelse(metadata.df$matatopy == 0, 0, 1))
log.cpm.extended[,"noMatatopy"] = as.factor(ifelse(metadata.df$matatopy == 1, 0, 1))
log.cpm.extended[,"boy"] = as.factor(ifelse(metadata.df$sex == 0, 0, 1))
log.cpm.extended[,"girl"] = as.factor(ifelse(metadata.df$sex == 1, 0, 1))
mod.matrix.extended = model.matrix(AD~., data = log.cpm.extended)

# --------------------------REPEATED CV-------------------------------------------------
#Lambda sequence for CV
lambda.seq = exp(seq(1,-5,length=200))
#Times to repeat the CV
n.repeat = 10
#Define vector of alphas to run cross validation
alphas =  seq(0.1, 1, 0.05)
#Define number of folds in CV
n.folds = 10



set.seed(345)
repeated.cv.results = repeat.cv.function(mod.matrix.extended, ad, alphas, lambda.seq, n.repeat, n.folds)
actual.rep.enet.model = glmnet(mod.matrix.extended[,-1], ad, family = "binomial",
                            alpha = repeated.cv.results$alpha, lambda = repeated.cv.results$lambda)


coeffs.actual.mod.rep = as.matrix(coef(actual.rep.enet.model))
names.final.coeffs.rep = str_remove_all(names(as.matrix(coeffs.actual.mod.rep)[as.vector(coeffs.actual.mod.rep) != 0,]), "`")
#xtable(as.matrix(coeffs.actual.mod.rep[as.vector(coeffs.actual.mod.rep) != 0,]))




# BOOTSTRAP REPEATED CV

#Define number of bootstrap samples
n.boot = 1000
set.seed(678)
boot.repeated.df = bootstrap.repeated.cv(log.cpm.extended, mod.matrix.extended, ad,
                                         alphas,lambda.seq, n.boot, n.repeat, n.folds)

#Save data from bootstrap
write.csv(boot.repeated.df[[2]], file = paste("Data/bootstrap-models-extended-repeated-folds-",
                                              n.folds,".csv",sep = ""), row.names = FALSE)
write.csv(boot.repeated.df[[1]], file = paste("Data/bootstrap-coefficients-extended-repeated-",
                                              n.folds, ".csv", sep = ""), row.names = FALSE)


#----------------------NESTED CV-------------------------------------------------




#Define number of folds for nested CV of lambda and alpha
n.folds.inner = 10
n.folds.outer = 10 
lambda.type = "lambda.1se"
alphas =  seq(0.1, 1, 0.05)



# RUN for actual model fit

set.seed(2345)

nested.cv.alpha.df = nested.cv.alpha(mod.matrix.extended, ad, n.folds.outer, n.folds.inner,
                                       alphas, lambda.type = lambda.type)
best.alpha.nested = nested.cv.alpha.df$alpha[nested.cv.alpha.df$deviance == min(nested.cv.alpha.df$deviance)]
final.model.cv = cv.glmnet(mod.matrix.extended[,-1], ad, family = "binomial", alpha = best.alpha.nested)


plot(final.model.cv$glmnet.fit, xvar = "lambda", label = TRUE)


coeffs.final.mod = as.matrix(coef(final.model.cv, s = lambda.type))
names.final.coeffs = str_remove_all(names(as.matrix(coeffs.final.mod)[as.vector(coeffs.final.mod) != 0,]), "`")
xtable(as.matrix(coeffs.final.mod[as.vector(coeffs.final.mod) != 0,], digits = 3))
#xtable(as.matrix(coeffs.actual.mod.rep[as.vector(coeffs.actual.mod.rep) != 0,]))



# BOOTSTRAP NESTED CV
set.seed(1235)
list.bootstrap = bootstrap.elasticnet(log.cpm.extended, mod.matrix.extended, response = ad, n.boot = 1000,
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
write.csv(list.bootstrap[[2]], file = paste("Data/bootstrap-models-extended-innerF",n.folds.inner,
                                            "-outerF", n.folds.outer,"-", lambda, "-new",
                                            ".csv",sep = ""), row.names = FALSE)
write.csv(list.bootstrap[[1]], file = paste("Data/bootstrap-coefficients-extended-innerF",
                                            n.folds.inner, "-outerF", n.folds.outer,"-",
                                            lambda,"-new", ".csv", sep = ""), row.names = FALSE)

