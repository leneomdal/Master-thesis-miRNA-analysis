source("Code//load_count_data.R")
library(edgeR)
library(DESeq2)

# NORMALIZE WITH EDGER TO CPM AND LOG CPM

#Make dge object
dge.full <- DGEList(counts=count.df)


# Filtering out miRNAs with cpm value less than 1 in 30 or less of the samples
keep = rowSums(cpm(dge.full)>1)>=10
dge.full$counts = dge.full$counts[keep,]


# Check for identical miRNAs that should be merged into one
#equal.mirnas = data.frame(V1 = c(0,0))
#i = 1
#for(m in 1:(nrow(dge.full$counts)-1)){
#  for(n in (m+1):nrow(dge.full$counts)){
#    if(m != n & all(dge.full$counts[m,] == dge.full$counts[n,])){
#      equal.mirnas[,i] = c(n,m)
#      i = i+1
#      print(paste(m,n))
#    }
#  }
#}

#Merge miRNA with identical expression and identical sequence
#Rename
rownames(dge.full$counts)[c(160, 426, 444, 449, 450, 461)] = c("miR-199a/b-3p", 
                                                        "miR-500a/b-5p", "miR-517a/b-3p", 
                                                        "miR-518d-5p-group", 
                                                        "miR-518e-5p-group", "miR-520b/c-3p")
#Delete duplicates
dge = dge.full[-c(161, 428, 445, 463, 470, 456, 458, 466, 467, 462),]

# normalize using TMM normalization
dge = calcNormFactors(dge, method = "TMM")
#Convert to cpm using TMM norm factors
cpm <- cpm(dge)


# log cpm
log.cpm = cpm(dge, log = TRUE)



# NORMALIZE WITH DESEQ USING REGULARIZED LOG
# (not used)

metadata.df$probiotic = as.factor(metadata.df$probiotic)
metadata.df$ad = as.factor(metadata.df$ad)
metadata.df$matatopy = as.factor(metadata.df$matatopy)

# Create DESeq matrix object
#ddsMat <- DESeqDataSetFromMatrix(countData = dge$counts,
#                                 colData = metadata.df,
#                                 design = ~ probiotic + ad)

#dds = ddsMat
# commented out cause uses filtering above
# Remove rows with zero or only 1 reads across all samples
#dds <- dds[ rowSums(counts(dds)) > 1, ]

# DESeq's regularized logtransformation
#rld <- rlog(dds, blind=TRUE)
# Plot first sample against second after log transformation
#plot(assay(rld)[,1:2], pch=16, cex=0.3)

# log2 transform after normalizing for comparison
#dds <- estimateSizeFactors(dds)
#plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1), pch=16, cex=0.3)



# Create groups
treatment = ifelse(metadata.df$probiotic ==1, "P", "nP")
outcome = ifelse(metadata.df$ad == 1, "AD", "nAD")
groups_p = as.factor(treatment)
groups_ad = as.factor(outcome)
groups = as.factor(paste0(treatment,"_", outcome))
groups_mad = as.factor(ifelse(metadata.df$matatopy == 1, "mAD", "no mAD"))



