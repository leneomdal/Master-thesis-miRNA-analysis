library(clValid)

intern.valid = clValid(log.cpm, 2:10, clMethods = c("hierarchical"),
                       validation = "internal", method = "ward", metric = "correlation")
?clValid
summary(intern.valid)
