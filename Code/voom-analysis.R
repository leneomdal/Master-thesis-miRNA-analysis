# USING VOOM limma TO ANALYSE MIRNA DATA

source("Code//normalize-and-filter.R")  # load data, perform filtering and 
                      #normalization using cpm limma and rlog DESeq
library(limma)
library(edgeR)

#FILTER USING EDGER
# dge <- DGEList(counts=count.df)
# 
# keep.exprs = filterByExpr(dge, group = groups_p)
# dge = dge[keep.exprs,, keep.lib.sizes = FALSE]
# dim(dge$counts)


# CONSIDER ONLY MIRNAS FROM ELASTICNET MODEL OF AD
#dge = dge$counts[names.final.coeffs[-1],]



# Define design matrices
design.matrix = model.matrix(~ probiotic + ad , data = metadata.df) # adjusting for ad
design.matrix.all = model.matrix(~probiotic + ad + matatopy, data = metadata.df)
design.matrix.intr = model.matrix(~probiotic + ad + probiotic*ad, data = metadata.df)



# VOOM method ( uses counts and not log cpm)
voom.weights = voomWithQualityWeights(dge, design.matrix.all, plot=FALSE)
voom.fit = lmFit(voom.weights, design.matrix.all)
eB.voom.fit = eBayes(voom.fit) 



#Results from voom fit
#topTable(eB.voom.fit, coef = 2, sort.by = "p")
#xtable(topTable(eB.voom.fit, coef = 2, sort.by = "p") %>% slice_head(n = 10), digits = 3)

#plotSA(eB.voom.fit, main="Final model voom: Mean-variance trend")


# RESULTING TOP MIRNA FROM VOOM
top.table = topTable(eB.voom.fit, coef = "probiotic1", sort.by = "p", number = 10)
topTable(eB.voom.fit, coef = 2, sort.by = "p", number = sum(top.table$P.Value<0.05))
topTable(eB.voom.fit, coef = 2, sort.by = "p")

top.mirna = topTable(eB.voom.fit, coef = 2, sort.by = "p", number = 1000)
top.mirna = rownames(top.mirna)








