library(readr)
library(dplyr)
#Read expression data
brca_nontumor <- read_delim("brca-rsem-fpkm-tcga.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
brca_tumor <- read_delim("brca-rsem-fpkm-tcga-t.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


#formatting
rn<-brca_nontumor$Hugo_Symbol
brca_nontumor<-brca_nontumor[,3:dim(brca_nontumor)[2]]
rownames(brca_nontumor)<-rn

rn<-brca_tumor$Hugo_Symbol
brca_tumor<-brca_tumor[,3:dim(brca_tumor)[2]]
rownames(brca_tumor)<-rn

#find correlation of tp53
tp53data<-brca_nontumor[which(rownames(brca_nontumor)=="TP53"),]

brca_nontumorCor<-cor(t(brca_nontumor))

x<-brca_nontumorCor[which(rownames(brca_nontumorCor)=="TP53"),]
sort(x, decreasing = T)
