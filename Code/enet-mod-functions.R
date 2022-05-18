

#                         ELASTICNET MODEL FUNCTIONS

library(glmnet)
library(stringr)

# --------------------Repeated CV------------------------------------------
repeat.cv.function = function(mod.matrix, response, alphas, lambda.seq, n.repeat, n.folds){
  rep.cv.df = data.frame(matrix(NA, nrow = length(alphas)*length(lambda.seq), ncol = 2 + n.repeat))
  colnames(rep.cv.df) = c("alpha", "lambda", rep("deviance", n.repeat))
  
  rep.alphas = c()
  for(i in seq_along(alphas)){
    rep.alphas = c(rep.alphas, rep(alphas[i], length(lambda.seq)))
  }
  rep.cv.df$alpha = rep.alphas
  rep.cv.df$lambda = rep(lambda.seq, length(alphas))
  
  #Define lambda folds, different for each repetition but equal for each alpha
  lambda.folds = list()
  seq.to.sample.folds = rep(1:n.folds, ceiling(nrow(mod.matrix)/n.folds))[1:nrow(mod.matrix)]
  for(r in seq_len(n.repeat)){
    lambda.folds[[r]] = sample(seq.to.sample.folds)
  }
  
  
  for(i in seq_len(n.repeat)){
    cv.mean.error = c()
    for(j in seq_along(alphas)){
      repeated.cv = cv.glmnet(mod.matrix, response, type.measure = "deviance", family = "binomial",
                              alpha = alphas[j], lambda = lambda.seq, foldid = lambda.folds[[i]])
      cv.mean.error = c(cv.mean.error, repeated.cv$cvm)
    }
    rep.cv.df[,i + 2] = cv.mean.error
    print(paste("repetition:", i))
  }
  
  # Find alpha lambda pair with smallest mean deviation over all repetitions
  rep.cv.df$mean.dev = rowMeans(rep.cv.df[,3:(ncol(rep.cv.df))])
  min.mean.dev = min(rep.cv.df$mean.dev)
  min.rep.cv = rep.cv.df[rep.cv.df$mean.dev == min.mean.dev,]
  
  # HVIS DET ER FLERE MED LIK? HVA GJØR MAN DA? TA DEN NEDERSTE?
  min.alpha = min.rep.cv$alpha[nrow(min.rep.cv)]
  min.lambda = min.rep.cv$lambda[nrow(min.rep.cv)]
  
  return(data.frame(alpha = min.alpha, lambda = min.lambda, mean.dev = min.mean.dev))
}

#---------------------------------------------------------------------------
#--------------------Bootstrap repeated CV----------------------------------
bootstrap.repeated.cv = function(log.cpm.ad, full.mod.matrix, response, alphas,lambda.seq, 
                                 n.boot, n.repeat, n.folds){
  boot.coeffs.df = data.frame(matrix(ncol = ncol(full.mod.matrix), nrow = 0))
  colnames(boot.coeffs.df) = colnames(full.mod.matrix)
  boot.enet.mods.df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(boot.enet.mods.df) = c("alpha", "lambda", "deviance")
  
  for(i in 1:n.boot){
    boot.index = sample(1:nrow(log.cpm.ad), size = nrow(log.cpm.ad), replace = TRUE)
    boot.data = log.cpm.ad[boot.index, ]
    boot.mod.matrix = model.matrix(AD~., data = boot.data)
    boot.ad = response[boot.index]
    
    boot.repeated.cv.results = repeat.cv.function(boot.mod.matrix, boot.ad, alphas, 
                                                  lambda.seq, n.repeat, n.folds)
    best.alpha = boot.repeated.cv.results$alpha
    best.lambda = boot.repeated.cv.results$lambda
    mean.dev = boot.repeated.cv.results$mean.dev
    
    boot.repeated.enet.model = glmnet(boot.mod.matrix[,-1], boot.ad, family = "binomial", 
                               alpha = best.alpha, 
                               lambda = best.lambda)
    
    
    #Save data from this bootstrap sample
    boot.coeffs = as.matrix(coef(boot.repeated.enet.model))
    boot.coeffs.df[i,] = boot.coeffs 
    boot.enet.mods.df[i,] = c(alpha = best.alpha, lambda = best.lambda, 
                              mean.deviance = mean.dev )
    print(paste("boostrap sample nr.", i))
  }
  return(list(boot.coeffs.df, boot.enet.mods.df))
}

#---------------------------------------------------------------------------

#---------------------Nested CV---------------------------------------------
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
    #print(paste("alpha:", alphas[j], "___________________________"))
    dev = 0
    for(i in seq_len(n.folds.outer)){                      # outer fold, nested CV 
      curr.test.fold = folds.outer == i
      curr.train.fold = folds.outer != i
      test.set = mod.matrix[curr.test.fold,]
      train.set = mod.matrix[curr.train.fold, ]
      
      cv.fit = cv.glmnet(train.set, response.var[curr.train.fold], type.measure = "deviance",
                         family = "binomial",
                         alpha = alphas[j], foldid = as.vector(list.foldid.inner[[i]]))              
      # inner fold, CV for lambda
      
      
      pred = predict(cv.fit, newx = test.set, s = lambda.type, 
                     type = "response")                    # predict on left out test set of CV
      
      dev = dev + find.dev(pred, as.numeric(as.character(response.var[curr.test.fold])))
      #print(paste("lambda:", signif(cv.fit$lambda.1se, digits = 3 ), "dev:", signif(find.dev(pred, as.numeric(as.character(response.var[curr.test.fold])))/sum(curr.test.fold), digits = 3), "pred:", pred[1], pred[2], pred[3], pred[4], "truth:", response.var[curr.test.fold][1], response.var[curr.test.fold][2], response.var[curr.test.fold][3], response.var[curr.test.fold][4] ))
    }
    
    cv.alpha.df[j,] = c(alphas[j], dev/nrow.mod.matrix)   # store average deviance across folds
    # together with its alpha value
  }
  toc()
  print("im here")
  return(cv.alpha.df)
}
#---------------------------------------------------------------------------
#--------------------Bootstrap nested CV------------------------------------

# Paired bootstrap NESTED CV
bootstrap.elasticnet = function(log.cpm.ad, full.mod.matrix, response, n.boot, n.folds.outer, 
                                n.folds.inner, alphas, lambda.type){
  boot.coeffs.df = data.frame(matrix(ncol = ncol(full.mod.matrix), nrow = 0))
  colnames(boot.coeffs.df) = colnames(full.mod.matrix)
  boot.enet.mods.df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(boot.enet.mods.df) = c("alpha", "lambda", "deviance")
  
  # define folds for CV of lambda after selecting alpha ( this should be equal for each bootstrap sample) why
  folds.lambda = sample(x = rep(1:n.folds.inner, ceiling(nrow(full.mod.matrix)/n.folds.inner)), 
                        size = nrow(full.mod.matrix), replace = FALSE)
  for(i in 1:n.boot){
    boot.index = sample(1:nrow(log.cpm.ad), size = nrow(log.cpm.ad), replace = TRUE)
    #boot.data = log.cpm.ad[boot.index, ]
    #boot.mod.matrix = model.matrix(AD~., data = boot.data)
    boot.mod.matrix = full.mod.matrix[boot.index,]
    boot.ad = response[boot.index]
    
    nested.cv.df = nested.cv.alpha(boot.mod.matrix, boot.ad, n.folds.outer, n.folds.inner, 
                                   alphas, lambda.type = lambda.type)
    
    #Find alpha with lowest mean deviance
    boot.alpha.best = nested.cv.df$alpha[nested.cv.df$deviance == min(nested.cv.df$deviance)]
    #If multiple alphas show equal deviance, choose the highest alpha
    if(length(boot.alpha.best) > 1){
      boot.alpha.best = boot.alpha.best[length(boot.alpha.best)]
    }
    if(boot.alpha.best >1){
      print("alpha is greater tha one here?")
    }
    
    boot.cv.enet = cv.glmnet(boot.mod.matrix[,-1], boot.ad, family = "binomial", alpha = boot.alpha.best, 
                             foldid = folds.lambda)
    
    if(lambda.type == "lambda.1se"){
      boot.lamba.best = boot.cv.enet$lambda.1se
      boot.dev = boot.cv.enet$cvm[boot.cv.enet$lambda == boot.cv.enet$lambda.1se]
    }else if(lambda.type == "lambda.min"){
      boot.lamba.best = boot.cv.enet$lambda.min
      boot.dev = boot.cv.enet$cvm[boot.cv.enet$lambda == boot.cv.enet$lambda.min]
    }else{
      print("Choose either lambda.min or lambda.1se for regularization")
    }
    
    #Save data from this bootstrap sample
    boot.coeffs = as.matrix(coef(boot.cv.enet, s = boot.lamba.best))
    boot.coeffs.df[i,] = boot.coeffs 
    boot.enet.mods.df[i,] = c(alpha = boot.alpha.best, lambda = boot.lamba.best, 
                              deviance = boot.dev)
    print(paste("boostrap sample nr.", i))
  }
  return(list(boot.coeffs.df, boot.enet.mods.df))
}

