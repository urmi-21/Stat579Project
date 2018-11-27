library(TCGAbiolinks)
library(dplyr)
library(DT)
library(plyr)
library(data.table)
packageVersion("TCGAbiolinks")
setwd("~/Downloads/RNAseqDB/RNAseqDB/data/tcgaMD")
#######################################################################################
#Function takes a df and expands it by unlisting elements at a column
expand<-function(df,colName){
  res<-data.frame()
  #for each row
  for(i in 1: dim(df)[1]){
    thisRow<-df[i, ! (colnames(df) %in% c(colName))]
    tempdf<-as.data.frame(df[i, c(colName)])
    #if list is empty skip that row
    if(dim(tempdf)[1]<1){
      next
    }
    #change colnames so they are unique
    colnames(tempdf)<-paste(paste(colName,".",sep = ""),colnames(tempdf),sep = "")
    #print(paste(i,colnames(tempdf)))
    newRow<-cbind(thisRow,tempdf)
    res<-bind_rows(res,newRow)
    
  }
  #print(res)
  return(res)
}

getjoinedBiospcClinc<-function(projName){
  print(paste("Downloading",projName))
  clinicalBRCA <- GDCquery_clinic(project = projName, type = "clinical")
  biospecimenBRCA <- GDCquery_clinic(project = projName, type = "Biospecimen")
  
  #rename all cols from clinical table with suffix clinical
  colnames(clinicalBRCA)<- paste0("clinical.",colnames(clinicalBRCA))
  
  #expand biospecimen data in the order portions, portions.analytes, portions.analytes.aliquots
  toUnpack<-c("portions", "portions.analytes", "portions.analytes.aliquots")
  for(s in toUnpack){
    biospecimenBRCA<-expand(biospecimenBRCA,s)
  }
  #add patient barcode to biospecimen data
  biospecimenBRCA<- biospecimenBRCA %>% mutate(clinical.bcr_patient_barcode=substr(submitter_id,1,nchar(as.character(submitter_id))-4))
  #join clinical and biospecimen
  finalJoined<-join(clinicalBRCA,biospecimenBRCA,by="clinical.bcr_patient_barcode")
  return(finalJoined)
}

##########################End Functions##########################################

#download and merge BRCA metadata
BRCAMetadata<-getjoinedBiospcClinc("TCGA-BRCA")
clinical <- GDCquery_clinic(project = "TCGA-UCS", type = "clinical")

tcgaProjList<-c("TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-UCEC","TCGA-UCS","TCGA-READ","TCGA-COAD","TCGA-LIHC","TCGA-HNSC","TCGA-ESCA","TCGA-PRAD","TCGA-STAD","TCGA-THCA","TCGA-LUAD","TCGA-LUSC","TCGA-KIRC","TCGA-KIRP","TCGA-KICH")
tcgaProjList<-c("TCGA-HNSC","TCGA-ESCA","TCGA-PRAD")

#mdList will have all data frames for rcgaProjList
mdListDF<-data.frame()
for(s in tcgaProjList){
  #mdList<-c(mdList,getjoinedBiospcClinc(s))
  if(dim(mdListDF)[1]<1){
    mdListDF<-getjoinedBiospcClinc(s)
  }else{
    print("joining")
    temp<-getjoinedBiospcClinc(s)
    mdListDF<-bind_rows(mdListDF,temp)  
  }
  
}

ulMD<-unlist(mdList)
mdJoined<-rbindlist(unlist(mdList))
n1<-colnames(t)
n2<-colnames(BRCAMetadata)
n3<-colnames(mdListDF)
setdiff(n2,n1)
#"updated_datetime" "submitter_id"     "created_datetime" "state" colnames are repeated


biospecimentest<- GDCquery_clinic(project = "TCGA-UCS", type = "Biospecimen")
clinicaltest<- GDCquery_clinic(project = "TCGA-UCS", type = "Clinical")

length(colnames(biospecimentest))
length(unique(colnames(biospecimentest)))

length(colnames(clinicaltest))
length(unique(colnames(clinicaltest)))

intersect(colnames(biospecimentest),colnames(clinicaltest))
