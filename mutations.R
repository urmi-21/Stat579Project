#process mutation data
library(dplyr)
library(DT)
library(plyr)
library(data.table)
library(maftools)
library("readr")
library(TCGAbiolinks)
packageVersion("TCGAbiolinks")

###########################################define functions##############################################################
#split maf file into categories by a given variable in the clinical metadata file
splitMafby<-function(clinicalData,by,mafData){
  #get tumor samps fo different values of by
  myList<-list()
  uniqVals<-clinicalData%>%select(by)%>%unique
  for(i in 1:nrow(uniqVals)){
    #i<-1
    s<-as.character(uniqVals[i,1])
    print(s)
    tList<-clinicalData %>% filter(clinicalData[,by]==s) %>% select(portions.analytes.aliquots.submitter_id)%>%unique
    thisData<-mafData%>%filter(Tumor_Sample_Barcode %in% tList$portions.analytes.aliquots.submitter_id)
    myList[[s]]<-thisData
  }
  
  return(myList)
}


##download gene mutation metadata
brcaMAF <- GDCquery_Maf("BRCA", pipelines = "varscan2")

#read clinical metadata
TCGAbrcaMetadata_reduced <- read_csv("TCGAbrcaMetadata_reduced.csv")

uv<-unique(TCGAbrcaMetadata_reduced$sample_type)

l1<-splitMafby(TCGAbrcaMetadata_reduced,"clinical.race",brcaMAF)

#for each item in list do calculations
lnames<-names(l1)
for(i in lnames){
  print((i))
  print(dim(l1[[i]]))
  #plot summary and save to pdf
  maf<-read.maf(l1[[i]],isTCGA = T)
  plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,showBarcodes=F)
}


clinicalData<-TCGAbrcaMetadata_reduced
by<-"clinical.race"
mafData<-brcaMAF

for(s in 1:nrow(uniqVals)){
  print(uniqVals[s,1])
  print("***")
}
