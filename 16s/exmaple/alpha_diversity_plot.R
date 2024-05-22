#alpha diversity function
alpha_plot <- function(alpha_diversity_qza,metadata,compaired,level,color_value){
  #Generates alpha_plot for the group column,if you want to view the diversity of other columns,please replace the other column's name to "group"
  
  p_list <- c("ggplot2", "ggsignif", "qiime2R","dplyr")
  
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
  data <- as.data.frame(read_qza(alpha_diversity_qza)$data)
  if (colnames(data)[1] == "shannon_entropy") {
    colnames(data)[1] <- "Shannon Index"
  } else if (colnames(data)[1] == "pielou_evenness") {
    colnames(data)[1] <- "Pielou's evenness"
  } else if (colnames(data)[1] == "faith_pd") {
    colnames(data)[1] <- "Faith's PD"
  } else if (colnames(data)[1] == "observed_features") {
    colnames(data)[1] <- "Observed features"
  }
  sample_columns <- grep("^sample", colnames(metadata), value = TRUE) #获取样品列
  data <- data[rownames(data) %in% metadata[,sample_columns],,drop = F] 
  name <- colnames(data)[1] 
  data$group <- metadata$group[match(rownames(data), metadata[,sample_columns])] #添加分组信息
  data$group <- factor(data$group, levels = level) #设置分组顺序
  #显著检验
  p <- ggplot(data, aes(group,data[,name])) +
    geom_boxplot(alpha = 0.7, 
                 size = 0.5,
                 width = 0.6) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin.., ), width = 0.1, color = "black") +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax.., ), width = 0.1, color = "black") +
    theme_classic() +
    geom_jitter( 
      size = 1.5, 
      alpha = 0.7, 
      aes(color = group),
      show.legend = FALSE,
      width = 0.02) +
    scale_color_manual(values = color_value) +
    xlab(NULL)+
    theme(axis.title = element_text(size = 8,face = "bold"),    
          axis.text = element_text(size = 8,face = "bold"),     
          legend.text = element_text(size = 2,face = "bold"),
          legend.title = element_text(size = 2,face = "bold")) +
    ylab(name) +
    geom_signif(comparisons = compaired,
                family =  "bold",
                size = 0.5,
                textsize= 2,
                step_increase = 0.05,
                map_signif_level = TRUE,
                test = t.test)
  return(p)
}

#exmaple
shannon <- "exmaple/shannon_vector.qza"
compaired <- list(c("subject-1", "subject-2"))
metadata <- read.table("exmaple/metadata.tsv",header = T,sep = "\t")
colnames(metadata)[7] <- "group"
p1 <- alpha_plot(shannon,metadata,compaired = compaired,level = c("subject-1","subject-2"),color_value = c('subject-1' = 'red', 'subject-2' = 'blue'))
ggsave("shannon.png",p1,width = 2, height = 2, units = "in")
