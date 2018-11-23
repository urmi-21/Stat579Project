library(TCGAbiolinks)
library(dplyr)
library(DT)
library(data.table)
library(plyr)
packageVersion("TCGAbiolinks")


query <- GDCquery(project = "TCGA-CHOL",  data.category = "Clinical", file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

queryB <- GDCquery(project = "TCGA-COAD",  data.category = "Biospecimen", file.type = "xml")
GDCdownload(queryB)

aliquot <- GDCprepare_clinic(queryB, clinical.info = c("aliquot"))
sample <- GDCprepare_clinic(queryB, clinical.info = c("sample"))
bio_patient <- GDCprepare_clinic(queryB, clinical.info = c("bio_patient"))
analyte <- GDCprepare_clinic(queryB, clinical.info = c("analyte"))
portion <- GDCprepare_clinic(queryB, clinical.info = c("portion"))
protocol <- GDCprepare_clinic(queryB, clinical.info = c("protocol"))
slide <- GDCprepare_clinic(queryB, clinical.info = c("slide"))

hbp<-as.data.frame(names(bio_patient))
hanalyte<-as.data.frame(names(analyte))
hportion<-as.data.frame(names(portion))
hprot<-as.data.frame(names(protocol))
hslid<-as.data.frame(names(slide))


join(aliquot,bio_patient)

length(unique(aliquot$bcr_patient_barcode))
length(unique(sample$bcr_patient_barcode))
length(unique(bio_patient$bcr_patient_barcode))
length(unique(analyte$bcr_patient_barcode))
length(unique(portion$bcr_patient_barcode))
length(unique(protocol$bcr_patient_barcode))
length(unique(slide$bcr_patient_barcode))
length(Reduce(intersect,list(aliquot$bcr_patient_barcode,sample$bcr_patient_barcode,bio_patient$bcr_patient_barcode,analyte$bcr_patient_barcode,portion$bcr_patient_barcode,protocol$bcr_patient_barcode,slide$bcr_patient_barcode)))

#remove duplicated rows from all tables
bio_patient_nr<-distinct(bio_patient)
aliquot_nr<-distinct(aliquot)
samp_nr<-distinct(sample)
analyte_nr<-distinct(analyte)
portion_nr<-distinct(portion)
protocol_nr<-distinct(protocol)
slide_nr<-distinct(slide)

#1 join aliquot with analyte, add analyte barcode in aliquot. for aliquot barcode TCGA-3L-AA1B-01A-01D-YYYY-23, analyte barcode is TCGA-3L-AA1B-01A-01D
#more info on barcodes https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
#tcga center codes https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/center-codes

aliquot_nr<-aliquot_nr%>% mutate(bcr_analyte_barcode=substr(bcr_aliquot_barcode,1,nchar(as.character(bcr_aliquot_barcode))-8))
j1<-join(analyte_nr,aliquot_nr,by="bcr_analyte_barcode")

#2 add bioportion barcode to j1 and join with bioportion
#portion: TCGA-3L-AA1B-01A-11 ; analyte: TCGA-3L-AA1B-01A-11
j1<-j1 %>% mutate(bcr_portion_barcode=substr(bcr_analyte_barcode,1,nchar(as.character(bcr_analyte_barcode))-1))
j2<-join(j1,portion_nr,by="bcr_portion_barcode")

#3 add biosample barcode to j2
#sample: TCGA-3L-AA1B-01A ; portion: TCGA-3L-AA1B-01A-11
j2<-j2 %>% mutate(bcr_sample_barcode=substr(bcr_portion_barcode,1,nchar(as.character(bcr_portion_barcode))-3))
j3<-join(j2,samp_nr,by="bcr_sample_barcode")

#finally join by biopatient
j4<-join(j3,bio_patient_nr,by="bcr_patient_barcode")


#download clinical data
queryC <- GDCquery(project = "TCGA-COAD",  data.category = "Clinical", file.type = "xml")
GDCdownload(queryC)
drug<-GDCprepare_clinic(queryC, clinical.info = "drug")
