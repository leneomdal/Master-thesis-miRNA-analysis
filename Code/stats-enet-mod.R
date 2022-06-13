

# Read saved data resulting from bootstrap nested cv of elasticnet model
#bootstrap.models.df = read.csv("Data//bootstrap-models-extended-innerF10-outerF10-lambda-1se-new.csv")
bootstrap.models.df = read.csv("Data//bootstrap-models-extended-repeated-folds-10.csv")
#bootstrap.models.df = read.csv("Data//bootstrap-models-repeated-folds-5.csv")
#View(bootstrap.models.df)
#bootstrap.coeffs.df = read.csv("Data//bootstrap-coefficients-extended-innerF10-outerF10-lambda-1se-new.csv")
bootstrap.coeffs.df = read.csv("Data//bootstrap-coefficients-extended-repeated-10.csv")
#bootstrap.coeffs.df = read.csv("Data//bootstrap-coefficients-repeated-5.csv")
colnames(bootstrap.coeffs.df)[1:(length(rownames(log.cpm))+1)] = c("(Intercept)", rownames(log.cpm))
#View(bootstrap.coeffs.df)



# extract coeffs included in at least one model

find.included.coeffs = function(bootstrap.coeffs.df){
  included.ind = c(NA)
  count = 1
  for(i in seq_len(ncol(bootstrap.coeffs.df))){
    if(sum(bootstrap.coeffs.df[,i] == 0) < 500){ # cut off for further plots
      included.ind[count] = i
      count = count +1
    }else{
    }
  }
  reduced.coeffs = bootstrap.coeffs.df[, included.ind]
  return(reduced.coeffs)
}


reduced.coeffs.df = find.included.coeffs(bootstrap.coeffs.df)
dim(reduced.coeffs.df)
dim(bootstrap.coeffs.df)


colnames(reduced.coeffs.df) = c(colnames(reduced.coeffs.df[1:(ncol(reduced.coeffs.df)-2)]), "matatopy", "noMatatopy")
reduced.coeffs.df = reduced.coeffs.df[,order(included.mirna.times.df$times, decreasing = TRUE)]
reduced.coeffs.df.long =  pivot_longer(reduced.coeffs.df[,-1], cols = colnames(reduced.coeffs.df[,-1]), 
                                       names_to = "miRNA", values_to = "coeffs")
reduced.coeffs.df.long$miRNA = factor(reduced.coeffs.df.long$miRNA, levels = colnames(reduced.coeffs.df), ordered = TRUE)


#final.coeffs.boot.long = pivot_longer(bootstrap.coeffs.df[,names.final.coeffs[-1]], 
#                                     cols = names.final.coeffs[-1], names_to = "miRNA", 
#                                       values_to = "coeffs")
# View(bootstrap.coeffs.df[,names.final.coeffs[-1]])





# Boxplot of bootstrap coeffs

ggplot(data = reduced.coeffs.df.long, aes(x = miRNA, y = coeffs)) + 
  geom_hline(yintercept = 0, color = "white") + geom_boxplot(outlier.size = 0.1, lwd = 0.001) + 
  theme( panel.grid = element_blank(),axis.text.x = element_text(angle = 90, 
                                                                       vjust = 0.5, hjust=1)) + 
  ylab("Coefficient estimate") + xlab("") +  ggtitle("Nested CV bootstrap samples")  +
  ylim(c(-4, 6))

# dat = data.frame(vec1 = c(2,2,2,2,1,1,1,1,-6,6,2,2,2,2,1,1,1,1,-6,6), vec2 = c(rep(1, 10), rep(2, 10)))
# ggplot(data = dat, aes(x = as.character(vec2), y = vec1)) + geom_boxplot()


#BARPLOT OF NUMBER OF TIMES EACH MIRNA IS INCLUDED

count.number.of.times.included = function(reduced.coeffs.df){
  included.times.df = data.frame(miRNA = names(reduced.coeffs.df),
                                 times = matrix(NA, nrow = ncol(reduced.coeffs.df), ncol = 1))
  count = 1
  for(i in 1:ncol(reduced.coeffs.df)){
    included.times.df$times[i] = sum(reduced.coeffs.df[,i] != 0)
    
  }
  return(included.times.df)
}
included.mirna.times.df = count.number.of.times.included(reduced.coeffs.df)
included.mirna.times.df$times = included.mirna.times.df$times
included.mirna.times.df$miRNA = factor(included.mirna.times.df$miRNA, levels = rev(included.mirna.times.df$miRNA), ordered = TRUE)
ggplot(data = included.mirna.times.df, aes(x = (1-times/1000), y = miRNA)) + 
         geom_bar(stat = "identity") +
  xlim(c(0,1)) + xlab("") + ggtitle("Proportion of boostrap estimated coefficients equal to 0") +
  theme(# remove the vertical grid lines
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank())


#histogram of alphas

ggplot(data = bootstrap.models.df, aes(x = alpha))  + geom_histogram(binwidth = 0.05, center = 0, color = "gray") +
  theme(panel.grid = element_blank()) + ggtitle("Alpha values chosen by repeated CV") + 
  xlab(expression(alpha)) +ylab(NULL) + ylim(c(0,450))

#historgram of lambdas

ggplot(data = bootstrap.models.df, aes(x = lambda))  + geom_histogram(binwidth = 0.04, center = 0, color = "gray") +
  theme(panel.grid = element_blank()) + ggtitle("Lambda values chosen by nested CV") + 
  xlab(expression(lambda)) +ylab(NULL) + ylim(c(0,380)) + xlim(c(-0.04,1.8))

