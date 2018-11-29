#process mutation data
library(dplyr)
library(DT)
library(plyr)
library(data.table)
library(maftools)
library("readr")
library(TCGAbiolinks)
library(ggplot2)
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

plotSummaryTofile<-function(mafList,fname){
  l1<-mafList
  #for each item in list do calculations
  lnames<-names(l1)
  plotList<-list()
  k<-0
  pdf(fname)
  for(i in lnames){
    print((i))
    print(dim(l1[[i]]))
    if(dim(l1[[i]])[1]<1){
      next
    }
    #plot summary and save to pdf
    maf<-read.maf(l1[[i]],isTCGA = T)
    plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = T, titvRaw = FALSE,showBarcodes=F)
    mtext(paste("Mutation summary by",i), outer=T,  cex=NA, line=-1.5,side = 3)
  }
  dev.off()
}


##############################################################################################################################

##download gene mutation metadata
brcaMAF <- GDCquery_Maf("BRCA", pipelines = "varscan2")

#read clinical metadata
TCGAbrcaMetadata_reduced <- read_csv("TCGAbrcaMetadata_reduced.csv")

tempList<-list()
tempList[["BRCADataset"]]<-brcaMAF
plotSummaryTofile(tempList,"mutationSummary.pdf")

l1<-splitMafby(TCGAbrcaMetadata_reduced,"clinical.race",brcaMAF)
plotSummaryTofile(l1,"mutationByRace.pdf")

l2<-splitMafby(TCGAbrcaMetadata_reduced,"clinical.primary_diagnosis",brcaMAF)
plotSummaryTofile(l2,"mutationByprimarydiagnosis.pdf")

l3<-splitMafby(TCGAbrcaMetadata_reduced,"clinical.gender",brcaMAF)
plotSummaryTofile(l3,"mutationBygender.pdf")

primaryDiagnosisFreq<-TCGAbrcaMetadata_reduced%>%select(clinical.primary_diagnosis)%>%group_by(clinical.primary_diagnosis)%>%count()%>%arrange(desc(freq))%>%mutate(logfreq=log10(freq))
primaryDiagnosisFreq$clinical.primary_diagnosis<-factor(primaryDiagnosisFreq$clinical.primary_diagnosis, levels = primaryDiagnosisFreq$clinical.primary_diagnosis)
ggplot(data=primaryDiagnosisFreq,aes(x=clinical.primary_diagnosis,y=logfreq,fill=freq))+geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

#group by diagnosis and gender
primaryDiagnosisFreq<-TCGAbrcaMetadata_reduced%>%select(clinical.gender,clinical.primary_diagnosis)%>%group_by(clinical.gender,clinical.primary_diagnosis)%>%count()%>%arrange(desc(freq))%>%mutate(logfreq=log10(freq))
ggplot(data=primaryDiagnosisFreq,aes(x=clinical.primary_diagnosis,y=logfreq,fill=clinical.gender))+geom_bar(stat = "identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#plot dnarna
ggplot(data=TCGAbrcaMetadata_reduced,aes(x=sample_type,fill=portions.analytes.analyte_type))+geom_bar()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

maf<-read.maf(l1[["white"]],isTCGA = T)
drugInteractions(maf = maf, fontSize = 0.75)
drugInteractions(genes = "DNMT3A", drugs = TRUE)

