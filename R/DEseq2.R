'''
DESeq2转录组差异表达分析
'''
setwd("~/PycharmProjects/dnad/GSE178341_Deseq2")
getwd()
mycounts<-read.csv("Epithelial cells_2.txt",row.names = 1, sep="\t",check.names = F)
#head(mycounts)
dim(mycounts)
mycounts_1<-mycounts[rowSums(mycounts) != 0,]
dim(mycounts_1)
mymeta<-read.csv("Epithelial cells_annotation_2.txt",stringsAsFactors = T,sep="\t")
#mymeta
#colnames(mycounts_1) == mymeta$sampleID


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=mycounts_1, 
                              colData=mymeta, 
                              design=~MSI)

#dds <- DESeq(dds)
dds <- DESeq(dds,sfType = "poscounts")
res <- results(dds)

head(res)
class(res)
res_1<-data.frame(res)
class(res_1)
head(res_1)
library(dplyr)
res_1 %>% 
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "UP",
    log2FoldChange <= -1 & padj <= 0.05 ~ "DOWN",
    TRUE ~ "NOT_CHANGE"
  )) -> res_2

table(res_2$group)

write.csv(res_2,file="result_fc2_1/Epithelial cells.csv",
          quote = F)

