source("Code//load_count_data.R")
library(edgeR)
library(DESeq2)

# NORMALIZE WITH EDGER TO CPM AND LOG CPM

#Make dge object
dge <- DGEList(counts=count.df)


#Filtering
keep = rowSums(dge$counts>3)>=3
dge$counts = dge$counts[keep,]
dim(dge)
sum(keep)


#keep = rowSums(dge$counts)>5
#sum(keep)
#dge$counts = dge$counts[keep,]
#dim(dge$counts)
#dim(count.df)

# normalize using TMM normalization
dge = calcNormFactors(dge, method = "TMM")
#Convert to cpm using TMM norm factors
cpm <- cpm(dge)
dim(cpm)


# log spm
log.cpm = cpm(dge, log = TRUE)



# NORMALIZE WITH DESEQ USING REGULARIZED LOG

metadata.df$probiotic = as.factor(metadata.df$probiotic)
metadata.df$ad = as.factor(metadata.df$ad)

# Create DESeq matrix object
ddsMat <- DESeqDataSetFromMatrix(countData = dge$counts,
                                 colData = metadata.df,
                                 design = ~ probiotic + ad)

dds = ddsMat
# commented out cause uses filtering above
# Remove rows with zero or only 1 reads across all samples
#dds <- dds[ rowSums(counts(dds)) > 1, ]
dim(dds)

# DESeq's regularized logtransformation
rld <- rlog(dds, blind=TRUE)
# Plot first sample against second after log transformation
#plot(assay(rld)[,1:2], pch=16, cex=0.3)

# log2 transform after normalizing for comparison
dds <- estimateSizeFactors(dds)
#plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1), pch=16, cex=0.3)



# Create groups
treatment = ifelse(metadata.df$probiotic ==1, "P", "nP")
outcome = ifelse(metadata.df$ad == 1, "AD", "nAD")
groups_p = as.factor(treatment)
groups_ad = as.factor(outcome)
groups = as.factor(paste0(treatment,"_", outcome))
groups_mad = as.factor(ifelse(metadata.df$matatopy == 1, "mAD", "no mAD"))

