# USING VOOM limma TO ANALYSE MIRNA DATA

source("Code//normalize-and-filter.R")  # load data, perform filtering and 
                      #normalization using cpm limma and rlog DESeq
library(limma)
library(edgeR)


# Define design matrices
design.matrix = model.matrix(~probiotic + ad , data = metadata.df) # adjusting for ad
design.matrix.all = model.matrix(~probiotic + ad + matatopy + sex + sib, data = metadata.df)
design.matrix.intr = model.matrix(~probiotic + ad + probiotic*ad, data = metadata.df)



# VOOM method ( uses counts and not log cpm)
voom.weights = voom(dge, design.matrix, plot=FALSE)
voom.fit = lmFit(voom.weights, design.matrix)
eB.voom.fit = eBayes(voom.fit) 


#Results from voom fit
#topTable(eB.voom.fit, coef = 2, sort.by = "p")
#xtable(topTable(eB.voom.fit, coef = 2, sort.by = "p") %>% slice_head(n = 10), digits = 3)

#plotSA(eB.voom.fit, main="Final model voom: Mean-variance trend")


# RESULTING TOP MIRNA FROM VOOM
top.table = topTable(eB.voom.fit, coef = 2, sort.by = "p", number = 10)
topTable(eB.voom.fit, coef = 2, sort.by = "p", number = sum(top.table$adj.P.Val<0.05))
top.mirna = topTable(eB.voom.fit, coef = 2, sort.by = "p", number = 20)
top.mirna = rownames(top.mirna)





