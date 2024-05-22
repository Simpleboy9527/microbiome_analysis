distance_qza <- "exmaple/bray_curtis_pcoa_results.qza"
metadata <- "exmaple/alpha_diversity/metadata.tsv"

beta_cpp <- function(distance_qza,asv_qza,metadata){
  p_list <- c("ggplot2", "qiime2R","dplyr","vegan")
  
  # Loop through each package in the list
  for (p in p_list) {
    # Check if the package is already installed
    if (!requireNamespace(p)) {
      # If not installed, install the package
      install.packages(p)
    }
    
    # Load the package into the current R session
    suppressWarnings(suppressMessages(library(p, character.only = TRUE)))
  }
 

  table <- as.data.frame(read_qza(asv_qza)$data)
  table <- as.data.frame(t(table))
  sample_columns <- grep("^sample", colnames(metadata), value = TRUE) #获取样品列
  #使table的行名按照metadata的样品列的顺序排列
  table$sampleid <- rownames(table)
  table <- table[match(metadata[,sample_columns],table$sampleid),]
  table <- table[,-ncol(table)]
  
  dune.div <- adonis2(table ~ group, data = metadata, permutations = 999, method="bray")
  
  dune.div
  
  dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)
  
  
  distance <- read_qza(distance_qza)$data

  x_label<-round(distance$ProportionExplained[1]*100,2) #获取PCoA1的解释度
  x_label
  y_label<-round(distance$ProportionExplained[2]*100,2) #获取PCoA2的解释度
  y_label
  distance_vectors <- distance$Vectors #获取PCoA坐标
  sample_columns <- grep("^sample", colnames(metadata), value = TRUE) #获取样品列
  distance_vectors <- distance_vectors[distance_vectors$SampleID %in% metadata[,sample_columns],,drop = F] 
  distance_vectors$group <- metadata$group[match(distance_vectors$SampleID, metadata[,sample_columns])] #添加分组信息
  distance_vectors
  
    p <- ggplot(distance_vectors, aes(x = round(distance_vectors$PC1,digits = 2), y = PC2)) +
      geom_point(size = 1.5, aes(color = group)) + 
      theme_classic()+# Use color aesthetic for points
      geom_vline(xintercept = 0, lty = "dashed") +
      geom_hline(yintercept = 0, lty = "dashed") +
      labs(x = paste0("PCoA1 ", x_label, "%"),
           y = paste0("PCoA2 ", y_label, "%")) +
      stat_ellipse(data = distance_vectors,
                   geom = "polygon",
                   aes(fill = group),
                   alpha = 0.3) +
      theme(axis.title = element_text(size = 8,face = "bold"),    
            axis.text = element_text(size = 8,face = "bold"),     
            legend.text = element_text(size = 8,face = "bold"),
            legend.title = element_text(size = 8,face = "bold"),
            panel.border = element_rect(fill=NA,color="black", size=.5, linetype="solid"))
  
 result <- list(p,dune.div)
}





