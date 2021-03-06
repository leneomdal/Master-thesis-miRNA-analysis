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
design.matrix = model.matrix(~ probiotic + ad + matatopy + sex, data = metadata.df) # adjusting for ad
design.matrix.all = model.matrix(~probiotic + ad + matatopy + ad*probiotic, data = metadata.df)
design.matrix.intr = model.matrix(~probiotic + ad + probiotic*ad, data = metadata.df)



# VOOM method ( uses counts and not log cpm)
voom.weights = voom(dge, design.matrix, plot=TRUE)
voom.fit = lmFit(voom.weights, design.matrix)
eB.voom.fit = eBayes(voom.fit, robust = FALSE) 

svg("C:\\Users\\Lene\\repositories\\master-thesis-latex\\figures\\mean-variance-trend-by-voom.svg")
voom.weights = voom(dge, design.matrix, plot=TRUE)
lines(c(-0.7,15.7), c(0.8017527,0.8017527), type="l", lty=2, col = "blue")
# Close the graphics device
dev.off() 



#Results from voom fit
#topTable(eB.voom.fit, coef = 2, sort.by = "p")
#xtable(topTable(eB.voom.fit, coef = 2, sort.by = "p") %>% slice_head(n = 10), digits = 3)

#plotSA(eB.voom.fit, main="Final model voom: Mean-variance trend")


# RESULTING TOP MIRNA FROM VOOM

top.table = topTable(eB.voom.fit, coef = "probiotic1", sort.by = "p", number = 50)
#xtable(top.table, digits = 3)
top.mirna = topTable(eB.voom.fit, coef = "probiotic1", sort.by = "p", number = 20)
top.mirna = rownames(top.mirna)



# vanilla limma
lcpm = cpm(dge, log = TRUE)
limma.fit = lmFit(lcpm, design.matrix)
fit.vanilla = eBayes(limma.fit, trend = FALSE)
#table for p vs nP
topTable(fit.vanilla)


s_0 = fit.vanilla$s2.prior
d_0 = fit.vanilla$df.prior



