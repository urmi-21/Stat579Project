library("biomaRt")
library(readr)
library(plyr)
listMarts()
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(ensembl)
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")


getSequence(id = c("673","7157","837"), 
            type="entrezgene",
            seqType="coding_gene_flank",
            upstream=100, 
            mart=ensembl) 



attList<-listAttributes(ensembl)

hsGeneData<-getBM(attributes = c("ensembl_gene_id","start_position","end_position","strand","transcript_count","percentage_gene_gc_content","gene_biotype","hgnc_symbol"),mart = ensembl)

#read tcga hgnc list and join with ensembl data
allCombined_20089_9149 <- read_csv("allCombined_20089_9149.csv")
head(colnames(allCombined_20089_9149))
colnames(hsGeneData)[which(colnames(hsGeneData)=="hgnc_symbol")]<-"Hugo_Symbol"

#join datasets
joinedDF<-join(hsGeneData,allCombined_20089_9149,by="Hugo_Symbol",type="full")
