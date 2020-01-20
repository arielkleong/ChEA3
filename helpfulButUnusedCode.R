## ENCODE. The col names without the GM
noGMEncode<-as.vector(gsub("_.*_","_",colnames(encode)))
encode3<-cbind(encode)
colnames(encode3)<-noGMEncode

## Find how many GMs were used
install.packages("gsubfn")
library(gsubfn)
middle<-as.vector(strapplyc(colnames(encode),"_(.*)_"))
GMs<-middle[grep('GM',middle)]
uniqueGMs<-unique(GMs)
length(unique(GMs))
length(unique(middle))

## Practice retrieving CREEDS IDs
sampleid<-character()
samplecount<-329
sampletf<-329
for(i in 1:2) {
  sampletf<-sampletf+1
  print(sampletf)
  
  r <- GET(paste0(CREEDS_URL, 'search'), query=list(q = tfs[i]))
  test <- fromJSON(httr::content(r, 'text'))
  if(length(test)!=0 && "hs_gene_symbol" %in% colnames(test)) {
    print("LOL")
    test = subset(test,test$hs_gene_symbol==tfs[i])
    test=subset(test,!(is.na(test["hs_gene_symbol"])))
    if(length(test$id)!=0) {
      sampleid<-append(sampleid, test$id)
      samplecount<-samplecount+1
    }
  }
  print("MEH")
}

## Make binary matrix, didn't reduce stuff to choose from so slower.
network<-matrix(0, length(genes),length(tfs))
network = as.data.frame(network)
colnames(network) = tfs
rownames(network) = genes
for(i in 1:length(tfs)) {
  chea_exp = chea3[,colnames(chea3)==tfs[i]]
  chea_exp = chea_exp[!is.na(chea_exp)]
  encode_exp = encode2[,colnames(encode2)==tfs[i]]
  encode_exp = encode_exp[!is.na(encode_exp)]
  for(j in 1:length(genes)){
    if(any(chea_exp == genes[j])&&any(encode_exp == genes[j])){
      network[j,i]=1
    }
    
  }
  
}

## Move to sources and targets table, first try
edges = data.frame("source","target")
colnames(edges) = c("x","y")
for(i in 1:length(colnames(network))){
  sources = rep(tfs[i],sum(network[,i]))
  targets = genes[network[,i]==1]
  temp_edges=data.frame(sources,targets)
  colnames(temp_edges)=c("x","y")
  edges = rbind(edges,temp_edges)
}

## Make network of tf-tf interactions. Did the same thing- called "reduced".
tfTargets<-subset(edges1, edges1[,2]%in%tfs)
tfTfInteractions<-write.table(tfTargets,file='tfTfInteractions.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

## To make reducedCREEDSCombined
count2<-0
for(i in 1) {
  count2<-count2+1
  print(count2)
  up<-stri_split_fixed(creeds[i,3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(creeds[i,4],",")
  down<-sapply(down,toupper)
  oneTF<-tfsInCREEDS[tfsInCREEDS[,1]==creeds[i,2],] # limit to the one tf of the row in creeds, with many targets
  if(length(oneTF[,1]) > 0) {
    upTF<-oneTF[oneTF[,2]%in%up,] # reduce to only targets in up regulated. maybe add info now, rbind.
    #print("problem? before")
    if(length(upTF[,1]) > 0) {
      upTF[,3]<-"up"
    }
    #print("problem after?")
    downTF<-oneTF[oneTF[,2]%in%down,]
    if(length(downTF[,1]) > 0) {
      downTF[,3]<-"down"
    }
    #print("here?")
    both<-rbind(upTF,downTF)
    #print("or here")
    if(length(both[,1]) > 0) {
      both[,4:6]<-NA
      both<-both[,c(1,3,2,4,5,6)]
      both[,4]<-creeds[i,1]
      both[,5]<-creeds[i,5]
      both[,6]<-creeds[i,6]

      colnames(both) = c("a","b","c","d","e","f")
      reducedCREEDSCombined<-rbind(reducedCREEDSCombined,both)
    }

  }
}

## Also for reducedCREEDSCombined
# Of the 75 transcription factors, only take the ones that are in CREEDS
reduced1<-reduced
reduced1[] <- lapply(reduced,as.character)
creeds1<-creeds

tfsInCREEDS<-reduced[reduced[,1]%in% creeds[,2],] #the tfs found in all 3 datasets
reducedCREEDSCombined<-data.frame("source", "interaction","target","series","cell type","perturbation", stringsAsFactors = FALSE)
colnames(reducedCREEDSCombined) = c("a","b","c","d","e","f")

count2<-0
for(i in 1:length(rownames(creeds))) {
  count2<-count2+1
  print(count2)
  up<-stri_split_fixed(creeds[i,3],",")
  up<-sapply(up,toupper)
  down<-stri_split_fixed(creeds[i,4],",")
  down<-sapply(down,toupper)
  oneTF<-tfsInCREEDS[tfsInCREEDS[,1]==creeds[i,2],] # limit to the one tf of the row in creeds, with many targets
  if(length(oneTF[,1]) > 0) {
    upTF<-oneTF[oneTF[,2]%in%up,] # reduce to only targets in up regulated. maybe add info now, rbind.
    #print("problem? before")
    if(length(upTF[,1]) > 0) {
      upTF[,3]<-"up"
    }
    #print("problem after?")
    downTF<-oneTF[oneTF[,2]%in%down,]
    if(length(downTF[,1]) > 0) {
      downTF[,3]<-"down"
    }
    #print("here?")
    both<-rbind(upTF,downTF)
    #print("or here")
    if(length(both[,1]) > 0) {
      both[,4:6]<-NA
      both<-both[,c(1,3,2,4,5,6)]
      both[,4]<-creeds[i,1]
      both[,5]<-creeds[i,5]
      both[,6]<-creeds[i,6]

      colnames(both) = c("a","b","c","d","e","f")
      reducedCREEDSCombined<-rbind(reducedCREEDSCombined,both)
    }

  }
}

## Test specific tf in CREEDS
response <- GET(paste0(CREEDS_URL, 'api'), query=list(id ='gene:1137'))
if (response$status_code == 200){
  response <- fromJSON(httr::content(response, 'text'))
  print(response)
}

## Labeling cell types
tFsCREEDSCellType[2:9,"g"]<-"hematopoietic cells"
tFsCREEDSCellType[10:11,"g"]<-"immune"
tFsCREEDSCellType[12:13,"g"]<-"liver"
tFsCREEDSCellType[14:18,"g"]<-"alveolar epithelial"
tFsCREEDSCellType[19,"g"]<-"immune"
tFsCREEDSCellType[20,"g"]<-"liver"
tFsCREEDSCellType[21:22,"g"]<-"immune"
tFsCREEDSCellType[23,"g"]<-"embryonic stem cells"
tFsCREEDSCellType[24,"g"]<-"neuron"
tFsCREEDSCellType[25:26,"g"]<-"immune"

## Getting a table for all tfs in CREEDS, and adding a column for human or mouse (organism)
install.packages("dplyr")
install.packages(c("nycflights13", "Lahman")) ## for examples

#tfsActuallyInCREEDS<-creeds[creeds[,2]%in% tfs,]
#tFsCREEDSCellType<-data.frame("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", stringsAsFactors = FALSE)

library(dplyr)
library(stringi)

tFsCREEDSCellType<-reducedCREEDSCombined
tFsCREEDSCellType[,c("g","h","i")]<-NA
tFsCREEDSCellType[1,7]<-"cell type"
tFsCREEDSCellType[1,8]<-"organism"
tFsCREEDSCellType[1,9]<-"dataset"
tFsCREEDSCellType[-1,9]<-"CREEDS"
tFsCREEDSCellType[1,4]<-"experiment (GEO series or PMID)"
tFsCREEDSCellType[1,5]<-"cell line"
tFsCREEDSCellType<-tFsCREEDSCellType[,c(1,2,3,4,5,7,6,8,9)]
#tFsCREEDSCellType[,7][tFsCREEDSCellType[,4]%in%reducedCREEDSCombinedHuman[,4]]<-"human" another way to do it
tFsCREEDSCellType<-transform(tFsCREEDSCellType, h = ifelse(tFsCREEDSCellType[,4]%in%reducedCREEDSCombinedHuman[,4],"human",h))

reducedCREEDSCombinedMouse<-anti_join(reducedCREEDSCombined,reducedCREEDSCombinedHuman)
tFsCREEDSCellType<-transform(tFsCREEDSCellType, h = ifelse(tFsCREEDSCellType[,4]%in%reducedCREEDSCombinedMouse[,4],"mouse",h))

## Get the cell lines I still need to identify
a<-creedsSep[!(creedsSep[,5]%in%allCellTypes[,1]),]

## Trying to process data so I can do it automatically
cellTypes2Revised<-read.table("C:/Users/maayanlab1/Downloads/cellTypes2 - cellTypes2.tsv",sep="\t",stringsAsFactors = FALSE)
uCT3Revised<-read.table("C:/Users/maayanlab1/Downloads/uniqueCellTypes3 - uniqueCellTypes3.tsv",sep="\t",stringsAsFactors = FALSE)

uCT3Revised[uCT3Revised[,1]=="Lin- Rosa26rtTA cells",2]<-"blood"

fill<-subset(uCT3Revised,is.na(uCT3Revised[,2])) # Make sure Lin-rosa had blood added; don't want to deal with spaces
fill1<-fill
fill[,1]<-sub(" (.*)","",fill[,1])

uCT3Revised<-uCT3Revised[!is.na(uCT3Revised[,2]),]

fillRevised<-read.table("C:/Users/maayanlab1/Downloads/fill - fill.tsv",sep="\t",stringsAsFactors = FALSE)
fillRevised[,1]<-fill1[,1]

allCellTypes<-rbind(cellTypes2Revised,uCT3Revised,fillRevised)
allCellTypes[719,1]<-creedsSepChEA[4596,5]

write.table(allCellTypes,"C:/Users/maayanlab1/Downloads/allCellTypes.tsv",sep="\t",row.names = FALSE,col.names = FALSE)

## Integrating cell types
cellTypes2<-read.csv("C:/Users/maayanlab1/Downloads/unique_cell_types - unique_cell_types.csv",stringsAsFactors = FALSE)
#cellTypes2<-as.data.frame(cellTypes2,stringsAsFactors=FALSE)
noNumbers<-as.vector(sub(".*? ","",cellTypes2[,1]))
cellTypes2[,1]<-noNumbers
for(i in 33:35) {
  cellTypes2[i,1]<-paste(cellTypes2[i,1],cellTypes2[i,3],sep=", ")
}
cellTypes2[,3]<-NULL
#cellTypes2<-cellTypes2[-234,] # it now skips row 234. 
colnames(cellTypes2)=c("cell lines","cell types")

cellTypes[,5]<-sub("TH1 ", "TH1",cellTypes[,5],fixed=TRUE) # space problem

for(i in 1:length(rownames(cellTypes2))) {
  cellTypes[which(cellTypes[,5]==cellTypes2[i,1]),6]<-as.character(cellTypes2[i,2])
  print(i)
}
NAs<-cellTypes[is.na(cellTypes[,6]),]

write.table(cellTypes2, file='C:/Users/maayanlab1/Downloads/cellTypes.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
write.table(cellTypes,file="C:/Users/maayanlab1/Downloads/cellTypesIn.csv",quote=FALSE,sep=',',row.names=FALSE,col.names=FALSE)

## Labeling GM cell types ENCODE
cellTypesNoGM<-cellTypes
cellTypesNoGM[,5][cellTypesNoGM[,5]=="GM18505"]<-"human"
cellTypesNoGM[,6][cellTypesNoGM[,5]=="GM18505"]<-"hematopoietic cells"
uniqueCellTypes<-as.data.frame(unique(cellTypes$`cell line`))

cellTypes<-read.table("C:/Users/maayanlab1/Downloads/cellTypesNoGMTSV.tsv", sep="\t",stringsAsFactors = FALSE)

## Parsing headers for ChEA
cellTypes<-tFsCREEDSCellType
# cellTypesCREEDSChEA<-data.frame("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", stringsAsFactors = FALSE)
# colnames(cellTypesCREEDSChEA)=c("a","b","c","d","e","f","g","h")

## Doesn't work as function for some reason, but separate seems to be okay.
 parseHeaders<- function(smaller, cellTypes, tfs) {
   for(i in 1:length(colnames(smaller))) { #length(colnames(smaller))
     header<-stri_split_fixed(colnames(smaller)[i],"_")
     #column<-smaller[!is.na(smaller[,i]) ,i]
     #column<-"bo"
     column<-smaller[smaller[,i]%in%tfs,i]
     if(length(column)>0) {
       headerRow<-matrix(NA,length(column),length(colnames(cellTypes)))
       headerRow<-as.data.frame(headerRow, stringsAsFactors = FALSE)
       colnames(headerRow)=c("a","b","c","d","e","f","g","h")
       headerNotOneItemList<-unlist(header)
       if(length(headerNotOneItemList)!=3) { ## specifically to distinguish between chea and encode
         headerRow[,1]<-headerNotOneItemList[1]
         headerRow[,3]<-as.character(column)
         headerRow[,4]<-headerNotOneItemList[2]
         headerRow[,5]<-paste(headerNotOneItemList[5:length(headerNotOneItemList)-1],sep=" ",collapse = " ")
         headerRow[,8]<-headerNotOneItemList[length(headerNotOneItemList)]
       }
       else {
         headerRow[,1]<-headerNotOneItemList[1]
         headerRow[,3]<-column
         headerRow[,5]<-headerNotOneItemList[2]
         headerRow[,8]<-headerNotOneItemList[3]
       }
       print(i)
       print(headerRow)
       cellTypes<-rbind(cellTypes,headerRow)
     }

   }
 }

 parseHeaders(chea2,cellTypes,tfs)
 parseHeaders(encode,cellTypes,tfs)

for(i in 1:length(colnames(chea2))) { #length(colnames(smaller))
  header<-stri_split_fixed(colnames(chea2)[i],"_")
  #column<-smaller[!is.na(smaller[,i]) ,i]
  #column<-"bo"
  column<-chea2[chea2[,i]%in%tfs,i]
  if(length(column)>0) {
    headerRow<-matrix(NA,length(column),length(colnames(cellTypes)))
    headerRow<-as.data.frame(headerRow, stringsAsFactors = FALSE)
    colnames(headerRow)=c("a","b","c","d","e","f","g","h","i")
    headerNotOneItemList<-unlist(header)
    
    headerRow[,1]<-headerNotOneItemList[1]
    headerRow[,3]<-as.character(column)
    headerRow[,4]<-headerNotOneItemList[2]
    headerRow[,5]<-paste(headerNotOneItemList[5:length(headerNotOneItemList)-1],sep=" ",collapse = " ")
    headerRow[,8]<-headerNotOneItemList[length(headerNotOneItemList)]
    headerRow[,9]<-"ChEA"
    
    print(i)
    print(headerRow)
    cellTypes<-rbind(cellTypes,headerRow)
  }
  
}
for(i in 1:length(colnames(encode))) { #length(colnames(smaller))
  header<-stri_split_fixed(colnames(encode)[i],"_")
  #column<-smaller[!is.na(smaller[,i]) ,i]
  #column<-"bo"
  column<-encode[encode[,i]%in%tfs,i]
  if(length(column)>0) {
    headerRow<-matrix(NA,length(column),length(colnames(cellTypes)))
    headerRow<-as.data.frame(headerRow, stringsAsFactors = FALSE)
    colnames(headerRow)=c("a","b","c","d","e","f","g","h","i")
    headerNotOneItemList<-unlist(header)
    
    headerRow[,1]<-headerNotOneItemList[1]
    headerRow[,3]<-column
    headerRow[,5]<-headerNotOneItemList[2]
    headerRow[,8]<-headerNotOneItemList[3]
    headerRow[,9]<-"ENCODE"
    
    print(i)
    print(headerRow)
    cellTypes<-rbind(cellTypes,headerRow)
  }
  
}

## Writing to file
 colnames(creedsSepChEA)=c("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", "source:target")
 colnames(creedsSepENCODE)=c("source", "interaction", "target", "experiment (GEO series or PMID)","cell line", "cell type", "perturbation", "organism", "source:target")
 write.table(creedsSepChEA,"C:/Users/maayanlab1/Downloads/creedsSepChEA.tsv",sep="\t",row.names = FALSE,col.names=TRUE)
 write.table(creedsSepENCODE,"C:/Users/maayanlab1/Downloads/creedsSepENCODE.tsv",sep="\t",row.names = FALSE,col.names=TRUE)
 
## Learning how to get cell lines from online
 #d<-read_html("http://web.expasy.org/cellosaurus/CVCL_0023")
 #tmp <- tempfile(fileext = ".xml")
 #write_xml(d, "C:/Users/maayanlab1/Downloads/alveolar.xml", options = "format")
 #readLines(tmp)
 #e<-htmlTreeParse("http://web.expasy.org/cellosaurus/CVCL_0023", useInternal=TRUE)
 
 e<-getURL("http://web.expasy.org/cellosaurus/CVCL_0023",ssl.verifypeer=FALSE)
 write(e,"C:/Users/maayanlab1/Downloads/alveolar.txt")
 
 ## Trying to actually use it:
 
 #fill[,2]<-as.character(fill[,2])
 
## Trying to figure out functions
 array<-matrix(0,2,2)
 double = function(number,array) {
   number<-number*2
   array[1,2]<-number
   return(array)
 }
 array<-double(2, array)
 
 ## Just seeing if I could've made yesterday easier
 
 cellLines_URL<-'https://scicrunch.org/resources/Cell%20Lines/search?q='
 e<-getURL("http://web.expasy.org/cellosaurus/CVCL_0023",ssl.verifypeer=FALSE)
 write(e,"C:/Users/maayanlab1/Downloads/alveolar.txt")
 
 fail<-getURL("https://scicrunch.org/resources/Cell%20Lines/search?q=pre-ipscs&l=pre-ipscs",ssl.verifypeer=FALSE)
 write(fail,"C:/Users/maayanlab1/Downloads/failedSearch.txt")
 
 
 #cellLine<-fill[i,1]
 # get_disease = function(cellLine, fill) {
 #   assign('fill',fill,envir=.GlobalEnv)
 #   cell_URL<- paste0('https://scicrunch.org/resources/Cell%20Lines/search?q=',cellLine,'&l=', cellLine, '&facet[]=Organism:Homo%20sapiens&facet[]=Organism:Mus%20musculus', collapse=NULL)
 #   b<-getURL(cell_URL,ssl.verifypeer=FALSE)
 #   write(b,"C:/Users/maayanlab1/Downloads/searchTest.txt")
 #   searchTestFile<-read_file("C:/Users/maayanlab1/Downloads/searchTest.txt")
 #   #print("a?")
 #   if(!grepl("We could not find",searchTestFile)) {
 #     #print("b?")
 #     disease<-strapplyc(b,"Disease:([^\"]*)")
 #     #print("c?")
 #     disease<-unlist(disease)
 #     #print(disease[1])
 #     #print(i)
 #     print(grep(cellLine, fill[,1]))
 #     fill[grep(cellLine, fill[,1]),2]<-disease[1]
 #     #print(fill[cellLine,2])
 #     #return(fill)
 #   }
 #   
 # }

## Trying to get URL
 #a <- GET(paste0(GEO_URL,'gse2527','\\)&retmax=10&usehistory=y', collapse=NULL))
 #b<-xmlParse(paste0(GEO_URL,'gse2527',')&retmax=10&usehistory=y', collapse=NULL))
 #b<-xmlParse(paste0(GEO_URL,'gse2527',')&retmax=10&usehistory=y', collapse=NULL))
 #b_data<-xmlToDataFrame(b)
 
## Notes:
 ## Get # tfs and targets in all three (function)
 # getCount = function(dataset) { same as dim()
 #   result<-c(length(colnames(dataset)),length(rownames(dataset)))
 #   return(result)
 # }
 
 # colnames and rownames (or dim) can give. 
 # chea tfs: 645, chea targets: 6208
 # encode tfs: 816, encode targets: 6921
 
 # They wanted all tf targets. 
 # creeds tfs: 575, creeds targets: 
 
 # get_organism = function(data_set_ids) {
 #   GEO_URL<-'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=('
 #   b<-getURL(paste0(GEO_URL,'gse2527',')&retmax=10&usehistory=y', collapse=NULL),ssl.verifypeer=FALSE)
 # }

 # for(i in 1:5) {
 #   fill<-get_disease(fill[i,1],fill)
 #   print(i)
 # }
 
 searchTest<-getURL(cell_URL,ssl.verifypeer=FALSE)
 write(searchTest,"C:/Users/maayanlab1/Downloads/searchTest.txt")
 searchTestFile<-read_file("C:/Users/maayanlab1/Downloads/searchTest.txt")
 
## Prepare to get whether mouse or human
 b<-getURL(paste0(GEO_URL,'gse2527',')&retmax=10&usehistory=y', collapse=NULL),ssl.verifypeer=FALSE)
 
 query_key<-strapplyc(b,"<QueryKey>(.*)</QueryKey>")
 web_env<-strapplyc(b,"<WebEnv>(.*)</WebEnv>")
 
 second_URL<-'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&query_key='
 
 c<-getURL(paste0(second_URL, query_key,'&WebEnv=',web_env, collapse=NULL),ssl.verifypeer=FALSE)
 organism<-strapplyc(c,"Organism:\t(.*?)\nType:")
 organism<-unlist(organism)
 organism<-organism[1]
 
## Fill in cell type for some of them
toFill<-subset(creedsSep,is.na(creedsSep[,6])) # found out there are 345 cell types to fill.

creedsSepENCODE<-creedsSep[creedsSep[,3]%in%encodeTfs,]
toFill2<-subset(creedsSepENCODE,is.na(creedsSepENCODE[,6]))
uniqueCellTypes3<-unique(c(toFill2[,5],toFill1[,5]))
uniqueCellTypes3<-as.data.frame(uniqueCellTypes3)

uniqueCellTypes3[,2]<-NA

## Matching cell lines to predict cell types

for(i in 1:length(rownames(uniqueCellTypes3))) {
  if(any(grepl(paste(cellTypes2[,1],collapse="|"),uniqueCellTypes3[i,1],ignore.case=TRUE))) {
    print(cat("i=",i))
    for(j in 1:length(rownames(cellTypes2))) {
      print(cat("j=",j))
      if(grepl(paste(cellTypes2[j,1],collapse="|"),uniqueCellTypes3[i,1],ignore.case=TRUE)) {
        uniqueCellTypes3[i,2]<-cellTypes2[j,2]
        break
      }
    }
    #print(grep(paste(pattern,collapse="|"),uniqueCellTypes2[i,1],ignore.case=TRUE))
    #uniqueCellTypes2[i,2]<-cellTypes2[grep(paste(pattern,collapse="|"),uniqueCellTypes2[i,1],ignore.case=TRUE),2]
  }
}

##  Writing to file reducedCREEDSCombined
reducedCREEDSCombinedTSV<-write.table(reducedCREEDSCombined,file='C:/Users/maayanlab1/Downloads/reducedCREEDSCombinedTSV.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

## Removing non-human ones: 
reducedCREEDSCombinedHuman<-reducedCREEDSCombined[reducedCREEDSCombined[,4]!="GSE24594",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE2527",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE2433",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE40273",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE6846",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE38375",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE21060",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE1566",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE46970",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE16974",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE32224",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE39009",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE10954",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE4356",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE27159",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE55272",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE39443",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE30323",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE12999",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE31354",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE47989",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE18383",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE33659",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE45941",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE5654",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE19923",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE11664",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE1948",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE34545",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE40296",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE23923",]
reducedCREEDSCombinedHuman<-reducedCREEDSCombinedHuman[reducedCREEDSCombinedHuman[,4]!="GSE28141",]

reducedCREEDSCombinedHumanTSV<-write.table(reducedCREEDSCombinedHuman,file='C:/Users/maayanlab1/Downloads/reducedCREEDSCombinedHumanTSV.tsv',quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)



for(i in 1:length(rownames(creedsSepChEA))) {
  print(i)
  creedsSepChEA[i,6]<-allCellTypes[which(creedsSepChEA[i,5]==allCellTypes[,1]),2]
}


for(i in 1:length(rownames(creedsSepENCODE))) {
  creedsSepENCODE[i,6]<-as.character(allCellTypes[which(creedsSepChEA[i,5]==allCellTypes[,1]),2])
}

creedsSepChEA<-transform(creedsSepChEA, h = ifelse(creedsSepChEA[,4]%in%reducedCREEDSCombinedHuman[,4],"human",h))
creedsSepChEA<-transform(creedsSepChEA, h = ifelse(creedsSepChEA[,4]%in%reducedCREEDSCombinedMouse[,4],"mouse",h))

creedsSepENCODE<-transform(creedsSepENCODE, h = ifelse(creedsSepENCODE[,4]%in%reducedCREEDSCombinedHuman[,4],"human",h))
creedsSepENCODE<-transform(creedsSepENCODE, h = ifelse(creedsSepENCODE[,4]%in%reducedCREEDSCombinedMouse[,4],"mouse",h))

## Getting unique series from creedsSepChEA and creedsSepENCODE
noOrganismChEA<-subset(creedsSepChEA,is.na(creedsSepChEA[,8]))
noOrganismENCODE<-subset(creedsSepENCODE,is.na(creedsSepENCODE[,8]))
noOrganism<-unique(rbind(noOrganismChEA,noOrganismENCODE))

## creedsSep after \t halted
for(i in 667516:length(rownames(creedsSep))) {
  print(i)
  creedsSep[i,6]<-allCellTypes[which(creedsSep[i,5]==allCellTypes[,1]),2] 
}

## grep finds every string that has the given as a substring
for(i in 440268:length(rownames(creedsSep))) {
  print(i)
  num<-grep(creedsSep[i,4],dataset_ids[,1],fixed=TRUE)
  if(length(num)) {
    creedsSep[i,8]<-dataset_ids[num,2]
  }
}