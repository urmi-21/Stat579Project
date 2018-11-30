library(readr)

#Read expression data
brca_nontumor <- read_delim("brca-rsem-fpkm-tcga.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
brca_tumor <- read_delim("brca-rsem-fpkm-tcga-t.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


#formatting
rn<-brca_nontumor$Hugo_Symbol
row.names(brca_nontumor)<-brca_nontumor$Hugo_Symbol
brca_nontumor<-brca_nontumor[,3:dim(brca_nontumor)[2]]

rownames(brca_tumor)<-brca_tumor$Hugo_Symbol
brca_tumor<-brca_tumor[,3:dim(brca_tumor)[2]]

#find correlation of tp53
which(rownames(brca_nontumor)=="tp53")
