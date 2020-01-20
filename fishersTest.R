install.packages("stats")

library(stats)
library(stringi)
library(dplyr)

## Find out how many rows the data frame should be
b<-apply(creeds,1,greatest)
max(b)

greatest<-function(x){
  up<-stri_split_fixed(x[3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(x[4],",")
  down<-sapply(down,toupper)
  
  return(length(up)+length(down))
}
## Format CREEDS to something that looks like ChEA and ENCODE

format<-function(d){
  up<-stri_split_fixed(d[3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(d[4],",")
  down<-sapply(down,toupper)
  name<-paste(d[2],d[1],d[5],d[6],collapse="_")
  temp<-data.frame(name=character(601),stringsAsFactors=FALSE)
  temp[,1]<-NA
  temp[1:length(up),1]<-up
  temp[(length(up)+1):(length(up)+length(down)),1]<-down
  colnames(temp)=name
  #print(typeof(temp))
  return(temp)
}

creedsFisher<-as.data.frame(do.call(cbind,(apply(creeds,1, function(d) format(d)))))
write.table(creedsFisher,"~/creedsFisher.tsv",col.names=T,sep="\t")

coexp<-read.table("~/coexp(1).tsv",header=T,stringsAsFactors = F)
coexp<-coexp500

## Do Fisher test
dataset1=chea1 
dataset2=coexp

if(!identical(dataset1,dataset2)) {
  colsCombo1<-combn(colnames(dataset1),2,simplify=FALSE)
  colsCombo2<-combn(colnames(dataset2),2,simplify=FALSE)
  
  colsCombo<-combn(union(colnames(dataset1),colnames(dataset2)),2,simplify=FALSE)
  colsCombo<-setdiff(colsCombo,colsCombo1)
  colsCombo<-setdiff(colsCombo,colsCombo2)
  
  genes<-unique(union(unlist(dataset1),unlist(dataset2)))%>%na.omit()
}

if(identical(dataset1,dataset2)) {
  colsCombo<-combn(colnames(dataset1),2,simplify=FALSE)
  genes<-unique(unlist(dataset1))%>%na.omit()
}

# Count time
ptm <- proc.time()

fishersJaccard<-function(d) { 
  
  tf1<-dataset1[[d[1]]]
  tf1<-tf1[!is.na(tf1)]
  
  tf2<-dataset2[[d[2]]]
  tf2<-tf2[!is.na(tf2)]
  
  jaccard<-(1-(length(intersect(tf1,tf2))/length(union(tf1,tf2))))
  
  contTable<-matrix(c((length(genes)-length(union(tf1,tf2))),length(setdiff(tf2,tf1)),length(setdiff(tf1,tf2)),length(intersect(tf1,tf2))),nrow=2)
  pVal<-fisher.test(contTable,alternative = "two.sided",conf.int=FALSE)
  
  row<-data.frame(d[1],d[2],pVal[[1]],jaccard,length(genes),length(tf1),length(tf2),length(intersect(tf1,tf2)),stringsAsFactors = FALSE)
  colnames(row)=c("tf1","tf2","pVal","jaccard","background genes","tf1 target genes","tf2 target genes","intersection size")
  return(row)
}

datasetPVal<-as.data.frame(do.call(rbind,(lapply(colsCombo, function(d) fishersJaccard(d))))) # PRACTICE

proc.time() - ptm

## Making the other columns
datasetPVal$match<-NA
datasetPVal$pValRank<-rank(datasetPVal$pVal,ties.method="max")
datasetPVal$jaccardRank<-rank(datasetPVal$jaccard,ties.method="max")