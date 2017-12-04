s1_LA <- read.delim("~/Desktop/coure project/DATA/GSM792454/s1_LA.txt")
s2_LA <- read.delim("~/Desktop/coure project/DATA/GSM792454/s2_LA.txt")
s3_LA <- read.delim("~/Desktop/coure project/DATA/GSM792454/s3_LA.txt")
s4_LA <- read.delim("~/Desktop/coure project/DATA/GSM792454/s4_LA.txt")                     
GSM54_LA<-as.matrix(cbind(s1_LA,s2_LA[,2],s3_LA[,2],s4_LA[,2]))
colnames(GSM54_LA)<-c("Gene","s1_LA","s2_LA","s3_LA","s4_LA")
s1_RA <- read.delim("~/Desktop/coure project/DATA/GSM792454/s1_RA.txt")
s2_RA <- read.delim("~/Desktop/coure project/DATA/GSM792454/s2_RA.txt")
s3_RA <- read.delim("~/Desktop/coure project/DATA/GSM792454/s3_RA.txt")
s4_RA <- read.delim("~/Desktop/coure project/DATA/GSM792454/s4_RA.txt") 
GSM54_RA<-cbind(s1_RA,s2_RA[,2],s3_RA[,2],s4_RA[,2])
colnames(GSM54_RA)<-c("Gene","s1_RA","s2_RA","s3_RA","s4_RA")
GSM54<-data.frame(GSM54_LA,GSM54_RA[,2:5])
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("EBSeq")
library(EBSeq)
data<-as.matrix(GSM54[,2:9])

###function ###normalization####
a<-cbind(median(s1_LA[,2]),median(s2_LA[,2]),median(s3_LA[,2]),median(s4_LA[,2]),median(s1_RA[,2]),median(s2_RA[,2]),median(s3_RA[,2]),median(s4_RA[,2]))
###transfer dataset to numeric data###
numfun<- function (x){
  for (i in 1:dim(x)[2]){
    dataNum[,i]<-c(as.numeric(x[,i]))
  }
  return (dataNum)
}
dataNum<-numfun(data)
row.names(dataNum)<-GSM54_LA[,1]

row.names(dataNum)<-GSM54_LA[,1]

####data normalization####
size<-MedianNorm(dataNum)
size<-as.matrix(size)


###apply ebseq function##
EBOut<-EBTest(Data=dataNum,Conditions=as.factor(rep(c("c1","c2"),each=4)),sizeFactors =size,maxround=20)
EBDEresult<-GetDEResults(EBOut,FDR=0.05)
head(EBDEresult$DEfound)
head(EBDEresult$PPMat)
str(EBDEresult$Status)
str(EBDEresult$DEfound)
###check for convergence###
EBOut$Alpha
EBOut$Beta
EBOut$P
which(is.element(rownames(EBDEresult$PPMat), EBDEresult$DEfound))
EBDEresult$PPMat[which(is.element(rownames(EBDEresult$PPMat), EBDEresult$DEfound)), ]
###Getting the result#####
###### Gene level DE analysis(two conditions)####

###calculate FC (Fold change of the raw data)###
GeneFC<-PostFC(EBOut)
PlotPostVsRawFC(EBOut,GeneFC)
par(mfrow=c(1,2))
QQP(EBOut)
DenNHist(EBOut)


