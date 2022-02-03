source("Code//normalize-and-filter.R")
library(pls)
library(plsRglm)


set.seed(5)

log.cpm.and.ad.df = as.data.frame((t(log.cpm)))
log.cpm.and.ad.df$ad = as.numeric(as.character(metadata.df$ad))
View(log.cpm.and.ad.df)
plsr.mod = plsr(ad~., data = log.cpm.and.ad.df, ncomp = 10, scale = FALSE, 
                center = FALSE,  validation = "CV", segments = 10, x = TRUE, 
                )
validationplot(plsr.mod, val.type = "MSEP")
?validationplot
?plsr
View(plsr.mod$x)

View(log.cpm)
View(count.df)
View(cpm)


?plsRglm

View(log.cpm.and.ad.df)
l = ncol(log.cpm.and.ad.df)
pls = plsRglm(ad~., data = log.cpm.and.ad.df, modele = "pls-glm-logistic", scaleX = FALSE, scaleY = FALSE, nt = 5)
pls$AIC
#cv.pls.mod = cv.plsRglm(ad~., data = log.cpm.and.ad.df[,(l-70):l],  modele = "pls-glm-logistic", scaleX = FALSE, scaleY = FALSE, NK = 1)

# Find first component
pls$Coeffs
pls$tt
pls$pp
pls$wwetoile
biplot(pls$tt, pls$pp[1:2,], col = c("blue", "red"))
?biplot

plot(pls$tt[,1], pls$tt[,2], col=ifelse(log.cpm.and.ad.df$ad, "blue", "red"))


pls$tt
pls$tt


sort(pls$pp[,1])
