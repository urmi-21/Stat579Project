library(readr)
library(dplyr)
library(ggplot2)
library(reshape)
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

#plot box plots of top mutated genes
topGenes<-c("PIK3CA","TP53","TTN","GATA3","CDH1","MAP3K1","MUC16","KMT2C","MUC4","PTEN")
nt<-as.data.frame(t(brca_nontumor[topGenes,]))
colnames(nt)<-topGenes

ggplot(data=nt)+ geom_boxplot(stat = "count")

ggplot(melt(nt), aes(x=factor(variable),y=value,fill=factor(variable)))+geom_boxplot()+scale_y_log10()+theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle=45,size = 15,face = "bold"),axis.text.y = element_text(size = 10,face = "bold"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), axis.title=element_text(size=12,face="bold"))+ scale_fill_discrete(name = "Gene")

