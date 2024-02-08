setwd("~/PycharmProjects/dnad/GSE139028_limma/")
options(stringsAsFactors = F)
rm(list=ls()) #清空变量
job <- "Macrophage" #设定项目名称

library(limma)
library(ggplot2) #用于绘制火山图
library(pheatmap) #用于绘制热图
library(dplyr)
library(ggVolcano)

### limma包差异表达分析 =================================================================
# 输入表达矩阵和分组文件 -------------------------------------------------------------
expr_data<-read.csv("Macrophage_limma.csv",header = T,row.names = 1,sep = ",",nrows = 7188) ##输入文件TPM原始值，行名是基因，列名是样本
expr_data <- expr_data[which(rowSums(expr_data)!=0),] #删除表达量为0的基因
expr_data = log2(expr_data + 1) #log化处理
expr_data[expr_data == -Inf] = 0 #将log化后的负无穷值替换为0
group<-read.csv("cells_info.csv",header = T,row.names = 1,sep = "\t") #输入文件，样本信息表，包含分组信息

# #构建分组矩阵--design ---------------------------------------------------------
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)

# #构建比较矩阵——contrast -------------------------------------------------------
contrast.matrix <- makeContrasts(burn-nor,levels = design) #根据实际的样本分组修改，这里对照组CK，处理组HT

# #线性拟合模型构建 ---------------------------------------------------------------
fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
#DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
#                       ifelse(DEG$logFC > 1, "up-regulated",
#                              ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
#write.table(table(DEG$regulate),file = paste0(job,"_","DEG_result_1_005.txt"),
#            sep = "\t",quote = F,row.names = T,col.names = T)
#write.table(data.frame(gene_symbol=rownames(DEG),DEG),file = paste0(job,"_","DEG_result.txt"),
#            sep = "\t",quote = F,row.names = F,col.names = T)

# 区分上下调基因 -----------------------------------------------------------------
#DE_1_0.05 <- DEG[DEG$P.Value<0.05&abs(DEG$logFC)>1,]
#upGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "up-regulated",]
#downGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "down-regulated",]
#write.csv(upGene_1_0.05,paste0(job,"_","upGene_1_005.csv"))
#write.csv(downGene_1_0.05,paste0(job,"_","downGene_1_005.csv"))

#tem1 <- head(rownames(upGene_1_0.05),1000) #以logFC差异倍数从大到小为序，提取前1000个基因名称
#tem2 <- as.data.frame(tem1) # 转化数据类型为数据框
#write.table(tem2$tem,file = paste0(job,"_","upgene_head100_name.txt"),
#            row.names = FALSE,col.names = FALSE)

#tem1 <- head(rownames(downGene_1_0.05),1000) #以logFC差异倍数从大到小为序，提取前1000个基因名称
#tem2 <- as.data.frame(tem1) # 转化数据类型为数据框
#write.table(tem2$tem,file = paste0(job,"_","downgene_head100_name.txt"),
#            row.names = FALSE,col.names = FALSE)

## 自定义筛选
foldChange = 4 # 自定义修改筛选参数
padj = 0.05 # 自定义修改筛选参数
All_diffSig <- DEG[(DEG$adj.P.Val < padj & (DEG$logFC > foldChange | DEG$logFC < (-foldChange))),]
#dim(All_diffSig)
write.csv(All_diffSig, paste0(job,"_","all_diffsig_filtered.csv"))  ##输出差异基因数据集

### 自定义筛选上调和下调的基因 ===================================================================
diffup <-  All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC > foldChange)),]
write.csv(diffup, paste0(job,"_","diffup_filtered.csv"))
diffdown <- All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC < -foldChange)),]
write.csv(diffdown, paste0(job,"_","diffdown_filtered.csv"))



# 火山图的绘制 ------------------------------------------------------------------
DEG$Genes <- rownames(DEG)
pdf(paste0(job,"_","volcano1.pdf"),width = 7,height = 7)
ggplot(DEG,aes(x=logFC,y=-log10(P.Value)))+ #x轴logFC,y轴adj.p.value
  geom_point(alpha=0.5,size=2,aes(color=regulate))+ #点的透明度，大小
  ylab("-log10(P.Value)")+ #y轴的说明
  scale_color_manual(values = c("blue", "grey", "red"))+ #点的颜色
  geom_vline(xintercept = c(-1,1),lty=4,col ="black",lwd=0.8)+ #logFC分界线
  geom_hline(yintercept=-log10(0.05),lty=4,col = "black",lwd=0.8)+ #adj.p.val分界线
  theme_bw()  #火山图绘制
dev.off()

# ggvolcano绘制另一种火山图 ----------------------------------------------------------

pdf(paste0(job,"_","volcano2.pdf"),width = 10,height = 10)
ggvolcano(data = DEG,x = "logFC",y = "P.Value",output = FALSE,label = "Genes",
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value",
          legend_position = "UR") #标签位置为up right
dev.off()
# 热图的绘制 -------------------------------------------------------------------
DEG_genes <- DEG[DEG$P.Value<0.05&abs(DEG$logFC)>1,]
DEG_gene_expr <- expr_data[rownames(DEG_genes),]
#DEG_gene_expr[is.infinite(DEG_gene_expr)] = 0
#DEG_gene_expr[DEG_gene_expr == -Inf] = 0
pdf(paste0(job,"_","pheatmap.pdf"))
pheatmap(DEG_gene_expr,
         color = colorRampPalette(c("blue","white","red"))(100), #颜色
         scale = "row", #归一化的方式
         border_color = NA, #线的颜色
         fontsize = 10, #文字大小
         show_rownames = F)
dev.off()

### 绘制火山图 ========================================================================
## 进行分类别
diffsig <- DEG
logFC <- diffsig$logFC
deg.padj <- diffsig$P.Value
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > 0.05 | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "Not"
data$group[(data$padj <= 0.05 & data$logFC > 1)] <-  "Up"
data$group[(data$padj <= 0.05 & data$logFC < -1)] <- "Down"
x_lim <- max(logFC,-logFC)

# 开始绘图
pdf(paste0(job,"_","volcano3.pdf"),width = 7,height = 6.5)  ## 输出文件
label = subset(diffsig,P.Value <0.05 & abs(logFC) > 0.5)
label1 = rownames(label)

colnames(diffsig)[1] = 'log2FC'
Significant=ifelse((diffsig$P.Value < 0.05 & abs(diffsig$log2FC)> 0.5), ifelse(diffsig$log2FC > 0.5,"Up","Down"), "Not")

ggplot(diffsig, aes(log2FC, -log10(P.Value)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  geom_text_repel(
    data = diffsig[diffsig$P.Value < 0.01&abs(diffsig$log2FC) > 6,],
    aes(label = Genes),
    size = 2,
    color = "black",
    segment.color = "black",
    show.legend = FALSE
  )+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  str(diffsig, max.level = c(-1, 1))+theme_bw()
dev.off()
