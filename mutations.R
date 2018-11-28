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
splitMafby<-function(clinicalData,by){
  
}


##download gene mutation metadata
brcaMAF <- GDCquery_Maf("BRCA", pipelines = "varscan2")

#read clinical metadata
