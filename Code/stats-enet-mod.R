

# Read saved data resulting from bootstrap nested cv of elasticnet model
bootstrap.models.df = read.csv("bootstrap-models.csv")
View(bootstrap.models.df)
bootstrap.coeffs.df = read.csv("bootstrap-coefficients.csv")
colnames(bootstrap.coeffs.df) = rownames(log.cpm)
View(bootstrap.coeffs.df)


# extract coeffs included in at least one model

find.included.coeffs = function(bootstrap.coeffs.df){
  included.ind = c(NA)
  count = 1
  for(i in seq_len(ncol(bootstrap.coeffs.df))){
    if(all(bootstrap.coeffs.df[,i] == rep(0, nrow(bootstrap.coeffs.df)))){}
    else{
      included.ind[count] = i
      count = count +1
    }
  }
  reduced.coeffs = bootstrap.coeffs.df[, included.ind]
  return(reduced.coeffs)
}


reduced.coeffs.df = find.included.coeffs(bootstrap.coeffs.df)

?pivot_longer
reduced.coeffs.df.long =  pivot_longer(reduced.coeffs.df, cols = colnames(reduced.coeffs.df), names_to = "miRNA", values_to = "coeffs")
View(reduced.coeffs.df.long)
