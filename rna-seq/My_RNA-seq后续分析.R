library(dplyr)       
library(ggplot2)     
library(pheatmap)
library(vegan)
library(DESeq2)
library(ggrepel)
library(ape)

#import feature-counts
feature_counts <- read.table('C:/Users/16236/Desktop/statistics/rna-seq/counts.txt',head=T,row.names=1)

my_count <- feature_counts[,-1:-5]

#filter 
result0 <- my_count[apply(my_count, MARGIN = 1, function(row) all(row > 0)),]
result <- my_count[apply(my_count, MARGIN = 1, function(row) all(row > 5)),]

#TPM标准化

kb <- feature_counts$Length/1000

RPK <- my_count/kb

sum <- colSums(RPK)/1000000

TPM <- matrix(0,40879,12)
for (i in 1:12) {
  TPM[,i] <- RPK[,i]/sum[i]
}
rownames(TPM) <- rownames(feature_counts)
colnames(TPM) <- c('D12','D16','D18','D19','D4','D7','J15','J16','J1','J4','J7','J8')


#DESeq2分析前准备，构建三个矩阵，表达矩阵，分组矩阵，差异比较矩阵
my_count <- as.matrix(my_count)

condition<-factor(c("D","D","D","D","D","D","J","J","J","J","J","J"))

coldata<-data.frame(row.names=colnames(my_count),condition) 

##构建dds矩阵 
dds<-DESeqDataSetFromMatrix(my_count,coldata,design=~condition)

head(dds) #查看构建好的矩阵

##进行差异分析
dds<-DESeq(dds) #对原始的dds进行标准化
vsd <- assay(vst(dds, blind=FALSE))#提取标准化结果
resultsNames(dds)  #查看结果名称
res<-results(dds) #用results函数提取结果，并赋值给res变量
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res, file="All_results.csv")


res_1 <- na.omit(res)
stat <- data.frame(rownames(res_1),res_1$stat,res_1$pvalue,res_1$padj)
stat <- stat[order(stat$res_1.stat),]
stat_filter <- stat[c(1:15, (nrow(stat) - 14):nrow(stat)), ]



# 打印结果
##热图

#
TPM_filter <- TPM[stat_filter$rownames.res_1.,]



p <- pheatmap(TPM_filter,scale = 'row',
              border="white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = T,
              angle_col = 45,
              show_rownames = T,
            
)

p



log10tpm <- log1p(TPM)

p <- pheatmap(t(log10tpm),
              cluster_cols = T, # 去掉横向、纵向聚类
              cluster_rows = T,
              show_rownames = T,
              show_colnames = F,
              border = T,
              border_color = "black")
p
  ##PCoA图




data <- t(TPM)

group <- data.frame(rownames(data),coldata[1])

distMatrix <- vegdist(data,method = "bray")

df.pcoa<-pcoa(distMatrix)

df.pcoa

df.plot<-data.frame(df.pcoa$vectors,group)

x_label<-round(df.pcoa$values$Relative_eig[1]*100,2)

y_label<-round(df.pcoa$values$Relative_eig[2]*100,2)

z_label<-round(df.pcoa$values$Relative_eig[3]*100,2)

a <- ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2,color=condition,shape=condition))+
  geom_point()+
  geom_text_repel(aes(label = rownames(df.plot)),size=3)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x=paste0("PCoA1 (",x_label,"%)"),
       y=paste0("PCoA2 (",y_label,"%)"))

b <- ggplot(data=df.plot,aes(x=Axis.1,y=Axis.3,color=condition,shape=condition))+
  geom_point()+
  geom_text_repel(aes(label = rownames(df.plot)),size=3)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x=paste0("PCoA1 (",x_label,"%)"),
       y=paste0("PCoA3 (",z_label,"%)"))

c <- ggplot(data=df.plot,aes(x=Axis.2,y=Axis.3,color=condition,shape=condition))+
  geom_point()+
  geom_text_repel(aes(label = rownames(df.plot)),size=3)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x=paste0("PCoA2 (",y_label,"%)"),
       y=paste0("PCoA3 (",z_label,"%)"))

p1 <- ggpubr::ggarrange(
  a,b,c,
  nrow = 2, ncol = 2,
  common.legend = TRUE, legend="right"
)
p1


#GO
library(clusterProfiler)
library(stringr)#基因ID转换
library(org.Mm.eg.db)
GO_database <- 'org.Mm.eg.db'

gene <- bitr(stat_filter$rownames.res_1.,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

#GO富集分析
enrich.go <- enrichGO(gene = gene$ENTREZID,  #待富集的基因列表
                      OrgDb = 'org.Mm.eg.db',  #指定物种的基因数据库，
                      keyType = 'ENTREZID',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                      qvalueCutoff = 1,  #指定 q 值阈值（可指定 1 以输出全部）
                      readable = FALSE)
write.table(enrich.go, 'enrich.go2.txt', sep = '\t', row.names = FALSE, quote = FALSE)
enrich.go <- as.data.frame(enrich.go)
MF <- enrich.go[enrich.go$ONTOLOGY == 'MF',]








#kegg 

gene_up_entrez <- as.character(na.omit(bitr(rownames(Sig_diff_gene[Sig_diff_gene$log2FoldChange > 1,]), #数据集
                                            fromType="SYMBOL", #输入格式
                                            toType="ENTREZID", # 转为ENTERZID格式
                                            OrgDb="org.Mm.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db" 

gene_entrez <- as.character(na.omit(bitr(rownames(Sig_diff_gene), #数据集
                                          fromType="SYMBOL", #输入格式
                                          toType="ENTREZID", # 转为ENTERZID格式
                                          OrgDb="org.Mm.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db" 

kegg_enrich <- enrichKEGG(gene  = gene_entrez,
                                  organism  = "mmu", #物种人类 hsa 小鼠mmu
                                  pvalueCutoff = 1,
                                  qvalueCutoff = 1)




kegg_enrich_results <- DOSE::setReadable(kegg_enrich, 
                                         OrgDb="org.Mm.eg.db", 
                                         keyType='ENTREZID')#ENTREZID to gene Symbol
export_data <-as.data.frame(kegg_enrich_results)
write.csv(as.data.frame(kegg_enrich_results),"kegg_reslut.csv")
