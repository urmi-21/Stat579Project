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

#download and combine maf files
getmaf<-function(projName){
  thisMaf<-GDCquery_Maf(projName, pipelines = "varscan2")
  return(thisMaf)
}

getmafs<-function(projList){
  allMaf<-data.frame()
  for(p in projList){
    temp<-GDCquery_Maf(p, pipelines = "varscan2")
    allMaf<-rbind(allMaf,temp)
  }
  return(allMaf)
}

############################################################################
tcgamafProjList<-c("BLCA","BRCA","CESC","UCEC","UCS","READ","COAD","LIHC","HNSC","ESCA","PRAD","STAD","THCA","LUAD","LUSC","KIRC","KIRP","KICH")
ucsmaf<-getmaf("UCS")
ucecmaf<-getmaf("UCEC")
uvmmaf<-getmaf("UVM")
mafs<-getmafs(c("UCS","UVM","BLCA"))

maf<-read.maf(mafs,isTCGA = T)
plotmafSummary(maf = maf)
