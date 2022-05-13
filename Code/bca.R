z0 = 1
a = 0.2

z.alpha = qnorm(0.025, lower.tail = TRUE)

compute.beta.quant = function(a, z0, z.alpha){
  betas = z0 + c( (z0 + z.alpha)/(1-a*(z0 + z.alpha)), (z0 - z.alpha)/(1-a*(z0 - z.alpha)))
  return(pnorm(betas))
}


#betas = b + c( (b + z.alpha)/(1-a*(b+z.alpha)), (b - z.alpha)/(1-a*(b - z.alpha)))


# Get bootrsapped results
bootstrap.models.df = read.csv("Data//bootstrap-models-repeated-folds-10.csv")
bootstrap.coeffs.df = read.csv("Data//bootstrap-coefficients-repeated-10.csv")
colnames(bootstrap.coeffs.df) = c("(Intercept)", rownames(log.cpm))


# define data frame to store values

df.confint.params = data.frame(matrix(nrow = ncol(bootstrap.coeffs.df), ncol = 6))
colnames(df.confint.params) = c("z0", "a", "beta1", "beta2","l.quantile", "u.quantile")


# Compute bias-correction factor z0
orig.coeffs = coeffs.actual.mod.rep
B = 1000
n.boot.coeffs.smaller = data.frame(n.smaller = NA)
for(i in 1:(length(orig.coeffs))){
  n = sum(bootstrap.coeffs.df[,i] < orig.coeffs[i])
  n.0 = sum(bootstrap.coeffs.df[,i] == orig.coeffs[i])
  n = n + n.0/2
  if(n == 0){
    n = 0.0001
  }
  n.boot.coeffs.smaller[i,1] = n
}
rownames(n.boot.coeffs.smaller) = colnames(bootstrap.coeffs.df)

z0 = qnorm(n.boot.coeffs.smaller$n.smaller/B)

df.confint.params$z0 = z0


# Estimating acceleration factor
n.samples = 60
calc.acc.factor = function(held.out.estimates,  n.samples){
  average.jackknife = sum(held.out.estimates)/n.samples
  a.top = 0
  a.bottom = 0
  for(i in 1:n.samples){
    a.top = a.top + (average.jackknife - held.out.estimates[i])^3
    a.bottom = a.bottom + (average.jackknife - held.out.estimates[i])^2
  }
  if(a.bottom == 0){
    a = 0
  } else{
    a = a.top/(6*a.bottom^(3/2))
  }
  return(a)
}


# # Estimate holding out one sample estimates for all coefficients
# df.held.out.estimates = data.frame(matrix( nrow = n.samples, ncol = length(orig.coeffs)))
# colnames(df.held.out.estimates) = colnames(bootstrap.coeffs.df)
# 
# 
# #Lambda sequence for CV
# lambda.seq = exp(seq(1,-5,length=200))
# #Times to repeat the CV
# n.repeat = 10
# #Define vector of alphas to run cross validation
# alphas =  seq(0.1, 1, 0.05)
# #Define number of folds in CV
# n.folds = 10
# 
# set.seed(54)
# for(i in 1:n.samples){
#   rep.cv.held.out = repeat.cv.function(mod.matrix[-i,], ad[-i], alphas, lambda.seq, n.repeat, n.folds)
#   rep.held.out.mod = glmnet(full.mod.matrix[-i,-1], ad[-i], family = "binomial",
#                          alpha = rep.cv.held.out$alpha, lambda = rep.cv.held.out$lambda)
#   coeffs.held.out = as.matrix(coef(rep.held.out.mod))
#   df.held.out.estimates[i,] = as.vector(coeffs.held.out)
# }
# write.csv(df.held.out.estimates, file = paste("Data/held-out-estimates.df", ".csv", sep = ""), row.names = FALSE)

df.held.out.estimates = read.csv("Data/held-out-estimates.df.csv", check.names = FALSE)

for(i in 1:(length(orig.coeffs))){
  a.curr =  calc.acc.factor(df.held.out.estimates[,i], n.samples)
  df.confint.params$a[i] = a.curr
  z0.curr = df.confint.params$z0[i]
  
  betas = compute.beta.quant(a.curr, z0.curr, z.alpha)
  beta1 = betas[1]
  beta2 = betas[2]
  df.confint.params$beta1[i] = beta1
  df.confint.params$beta2[i] = beta2
  df.confint.params$l.quantile[i] = quantile(bootstrap.coeffs.df[,i], probs = beta1)
  df.confint.params$u.quantile[i] = quantile(bootstrap.coeffs.df[,i], probs = beta2)
}

write.csv(df.confint.params, file = paste("Data/confint-boot-params", ".csv", sep = ""), row.names = FALSE)

# Check for miRNAs with confidence inetrvals not containing 0
for(i in 1:length(orig.coeffs)){
  if(df.confint.params$l.quantile[i] > 0 & df.confint.params$u.quantile[i] > 0){
    print(rownames(df.confint.params)[i])
  }
  if(df.confint.params$l.quantile[i] < 0 & df.confint.params$u.quantile[i] < 0){
    print(rownames(df.confint.params)[i])
  }
}

rownames(df.confint.params) = colnames(bootstrap.coeffs.df)

df.confint.params.included = df.confint.params[names.final.coeffs,]


# plot confidence interval

ggplot(data = final.coeffs.boot.long[], aes(x = miRNA, y = coeffs)) + geom_boxplot()

