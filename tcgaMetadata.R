library(TCGAbiolinks)
library(dplyr)
library(DT)
library(plyr)
packageVersion("TCGAbiolinks")

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
    #for(j in 1: dim(tempdf)[1]){
    #convert to dataframe in case there is only single column
    #res<-bind_rows(res,newRow)
    #}
  }
  #print(res)
  return(res)
}

getjoinedBiospcClinc<-function(projName){
  clinicalBRCA <- GDCquery_clinic(project = projName, type = "clinical")
  biospecimenBRCA <- GDCquery_clinic(project = projName, type = "Biospecimen")
  #expand biospecimen data in the order portions, portions.analytes, portions.analytes.aliquots
  toUnpack<-c("portions", "portions.analytes", "portions.analytes.aliquots")
  for(s in toUnpack){
    biospecimenBRCA<-expand(biospecimenBRCA,s)
  }
  #add patient barcode to biospecimen data
  biospecimenBRCA<- biospecimenBRCA %>% mutate(bcr_patient_barcode=substr(submitter_id,1,nchar(as.character(submitter_id))-4))
  #join clinical and biospecimen
  finalJoined<-join(clinicalBRCA,biospecimenBRCA,by="bcr_patient_barcode")
  return(finalJoined)
}

##########################End Functions##########################################

tcgaProjList<-c("")
BRCAMetadata<-getjoinedBiospcClinc("TCGA-BRCA")








