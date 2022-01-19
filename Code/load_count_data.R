library(limma)
library(tidyverse)

path = "C:\\Users\\Lene\\Documents\\Skole\\Prosjektoppgave\\project-thesis-mirna\\new-data\\"

#Read the data
count.data <- read.csv(paste0(path, "MatureMatrix.csv"),sep="\t", header =T,row.names=1)
sample.sheet<-read.table(paste0(path,"SampleSheet.txt"),sep="\t", header=T)

#Remove "X" from colnames
colnames(count.data) = str_remove(colnames(count.data), "[X]")

#Check for missing values
sum(is.na(sample.sheet) == TRUE)# Two samples missing value for sens2yrs, row 16 and 64

#Check for same samples and same ordering in count matrix and sample sheet
count.data<-count.data[,order(as.numeric(colnames(count.data)))]
sum(colnames(count.data) == sample.sheet$samples)



#Find only those included in miRNA sequencing project for 10 day samples, mirna = 1
ten.days.sample.included.nr = sample.sheet$bm_sample_no1[sample.sheet$day10 == 1]
n.samples.10.days = length(ten.days.sample.included.nr)


# Check that all samples are present
#for (i in seq(from = 1,to = length(ten.days.sample.included.nr))) {
#  if(toString(ten.days.sample.included.nr[i] %in% colnames(count.data))){
#    print("im present")
#  }
#  else {print(ten.days.sample.included.nr[i])}
#}



#Define data frame containing only the samples we want
ten.days.sample.df = count.data[,which(sample.sheet$day10 == 1)]


#Check for missing values in these extracted samples
apply(is.na(ten.days.sample.df), 2, which)
sum(is.na(ten.days.sample.df))


#Data frame only for needed metadata
ten.days.meta.data = sample.sheet[sample.sheet$day10 ==1,]
#nrow(ten.days.meta.data)

for(i in 1:ncol(ten.days.sample.df)){
  if(names(ten.days.sample.df)[i] != ten.days.meta.data$samples[i]){
    print("not present!")
  }
}

count.df = ten.days.sample.df
metadata.df = ten.days.meta.data
#View(metadata.df)
#ncol(count.df)
#nrow(metadata.df)


#Change col name of atopic dermatitis collumn
colnames(metadata.df)[colnames(metadata.df) == 'ad_ukwp2yr'] = 'ad'
