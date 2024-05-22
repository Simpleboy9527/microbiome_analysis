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

colnames(my_count) <- c('D12','D16','D18','D19','D4','D7','J15','J16','J1','J4','J7','J8')

#TPM标准化

kb <- feature_counts$Length/1000  #基因长度单位kb

RPK <- my_count/kb

sum <- colSums(RPK)/1000000

TPM <- matrix(0,40879,12)
for (i in 1:12) {
  TPM[,i] <- RPK[,i]/sum[i]
}
rownames(TPM) <- rownames(feature_counts)
colnames(TPM) <- c('D12','D16','D18','D19','D4','D7','J15','J16','J1','J4','J7','J8')

TPM_filtered  <- TPM[rowSums(TPM > 15, na.rm = TRUE) == ncol(TPM), ]

TPM_log10 <- log10(TPM_filtered+1) 


annotation_col = data.frame(Group = factor(rep(c("Control", "Nisin"), each = 6)))
rownames(annotation_col) = colnames(TPM)
pheatmap(TPM_log10,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_cols = T,
         cluster_rows = TRUE,
         cellwidth = 25,
         cellheight = 0.07,
         filename = "基因家族热图.pdf",
)

#DESeq2分析前准备，构建三个矩阵，表达矩阵，分组矩阵，差异比较矩阵
my_count1 <- as.matrix(my_count)

my_count1 <- my_count1[rownames(my_count1) %in% rownames(TPM_filtered),]

condition<-factor(c("D","D","D","D","D","D","J","J","J","J","J","J"))

coldata<-data.frame(row.names=colnames(my_count1),condition) 

##构建dds矩阵 
dds<-DESeqDataSetFromMatrix(my_count1,coldata,design=~condition)

head(dds) #查看构建好的矩阵

##进行差异分析
dds<-DESeq(dds) #对原始的dds进行标准化
vsd <- assay(vst(dds, blind=FALSE))#提取标准化结果
resultsNames(dds)  #查看结果名称
res<-results(dds) #用results函数提取结果，并赋值给res变量
res = res[order(res$pvalue),]
head(res)
summary(res)

#根据stat统计量(stat,即wald test)进行筛选上下调基因
res_1 <- na.omit(res)
stat <- data.frame(rownames(res_1),res_1$stat,res_1$pvalue,res_1$padj)
stat <- stat[order(stat$res_1.stat),]
stat_filter <- stat[c(1:15, (nrow(stat) - 14):nrow(stat)), ]


#热图
TPM_filter <- TPM[stat_filter$rownames.res_1.,]

TPM_filter <- TPM_filter[,c("D19","D18","D16","D4","D12","D7","J16","J15","J7","J8","J1","J4")]

p <- pheatmap(TPM_filter,scale = "row",
              border="white",  # 设置边框为白色
              cluster_cols = F, #去掉横向聚类、纵向聚类
              cluster_rows = F, 
              angle_col = 45,  #横坐标倾斜45°
              show_rownames = T, #显示gene名
              
)

p


#火山图
data <- as.data.frame(na.omit(res))

data %>% mutate(regulate = case_when(log2FoldChange>1&padj<0.05 ~ "up",
                                    log2FoldChange<(-1)&padj<0.05 ~ "down",
                                    TRUE ~ "ns")) ->data


p1 <- ggplot(data,aes(log2FoldChange,-log10(padj),color = regulate))+
  geom_point()+
  scale_color_manual(values = c("blue", "gray", "red"))+
  geom_hline(yintercept =  -log10(0.05),colour= "red", linetype = "dashed")
p1
ggsave("huoshantu.pdf",p1)

p2<- ggplot(data,aes(log2FoldChange,-log10(pvalue),color = regulate))+
  geom_point()+
  scale_color_manual(values = c("blue", "gray", "red"))+
  geom_hline(yintercept =  -log10(0.05),colour= "red", linetype = "dashed")
p2



#组间

#样本相似性
dat_500 <- TPM_filtered[names(sort(apply(TPM_log10,1,mad),decreasing = T)[1:500]),]#取高表达量前500基因
M <- cor(dat_500,method = "spearman")

p3 <-pheatmap::pheatmap(M,
                        show_rownames = T,
                        angle_col=45,
                        annotation_col = annotation_col) 
p3
ggsave(p3,filename = 'check_cor_top500.pdf',width = 7.5,height =6)


p4 <- pheatmap::pheatmap(cor(my_count,method = "spearman"),
                         show_rownames = T,
                         angle_col=45,
                         annotation_col = annotation_col) 
p4
ggsave(p4,filename = 'check_cor.pdf',width = 7.5,height =6)

p5 <- pheatmap::pheatmap(cor(TPM_filtered,method = "spearman"),
                         show_rownames = T,
                         angle_col=45,
                         annotation_col = annotation_col) 
p5
ggsave(p5,filename = 'check_cor_filtered.pdf',width = 7.5,height =6)


#组内多样性
groupa <- TPM_log10[,1:6]
groupb <- TPM_log10[,7:12]

t_test_result <- apply(groupa, 1, function(row) {
  result <- t.test(row)
  return(c(t_statistic = result$statistic, p_value = result$p.value))
})

# 输出每行的 t 统计量和 p 值
print(t_test_result)


TPM_log10

