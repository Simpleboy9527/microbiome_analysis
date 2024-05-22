taxa_barplot_cpp <- function(metadata, SVs, taxonomy, name, order) {
  # List of packages to check and install
  p_list <- c("ggplot2", "tidyverse", "qiime2R","dplyr","ggsci")
  
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
  
  
  # Read feature table (SVs) and taxonomy data
  SVs_data <- as.data.frame(read_qza(SVs)$data)
  taxonomy_data <- read_qza(taxonomy)$data %>% parse_taxonomy()
  
  # Summarize taxa based on specified name
  taxasums <- summarize_taxa(SVs_data, taxonomy_data)[[name]]
  i <- match(name, colnames(taxonomy_data))
  
  # Rename rows based on taxonomic hierarchy
  new_row_names <- sapply(strsplit(row.names(taxasums), ";"), function(x) paste(x[i], collapse = ";"))
  taxasums$taxa <- new_row_names
  taxasums$taxa[taxasums$taxa == " NA"] <- " Others"
  
  # Order taxa by abundance
  taxasums <- taxasums[order(rowSums(taxasums[, 1:(ncol(taxasums)-1)]), decreasing = TRUE), , drop = TRUE]
  taxasums <- rbind(taxasums[taxasums$taxa != " Others", ], taxasums[taxasums$taxa == " Others", ])
  
  # If more than 15 taxa, group less abundant ones as "Others"
  if (nrow(taxasums) > 15) {
    taxasums$taxa[16:nrow(taxasums)] <- " Others"
  }
  
  # Create a data frame with aggregated taxa abundances
  result_df <- taxasums %>%
    group_by(taxa) %>%
    summarise(across(colnames(taxasums)[-ncol(taxasums)], sum))
  
  result_df <- as.data.frame(result_df)
  row.names(result_df) <- result_df$taxa
  result_df <- result_df[, -1, drop = FALSE]
  
  # Check if data needs normalization
  column_sums_result_df <- colSums(result_df)
  are_column_sums_equal <- all(column_sums_result_df == column_sums_result_df[1])
  
  if (are_column_sums_equal) {
    # Normalize data if column sums are equal
    result_df_prop <- result_df / column_sums_result_df[1]
    
    # Create a wide matrix for plotting
    wide_matrix <- rownames_to_column(result_df_prop, var = "level")
    
    # Convert wide matrix to long matrix
    long_matrix <- gather(wide_matrix, key = "Sample", value = "Value", -level)
    
    # Match group information from metadata
    long_matrix$group <- metadata$group[match(long_matrix$Sample, metadata$sample.id)]
    
    # Set color palettes
    col = pal_d3("category20")(20)
    col2 = pal_d3("category20", alpha = 0.5)(20)
    mypal = c(col, col2[-8])
    
    # Order taxa and groups for plotting
    order_x = apply(result_df_prop[, 1:ncol(result_df_prop)], 1, sum)
    order_x = order_x[order(order_x, decreasing = TRUE)]
    site <- match(" Others", names(order_x))
    
    if (!is.na(site) && site != 1) {
      order_x <- order_x[c(site, 1:(site-1), (site+1):length(order_x))]
    }
    
    long_matrix$level = factor(long_matrix$level, levels = names(order_x), ordered = TRUE)
    long_matrix$group = factor(long_matrix$group, levels = order, ordered = TRUE)
    
    # Create and customize the bar plot using ggplot2
    temp_plot <- ggplot(long_matrix, aes(x = Sample, weight = Value, fill = level)) +
      geom_bar(aes(x = factor(Sample, levels = unique(long_matrix$Sample))), position = "stack") +
      theme_classic() +
      scale_fill_manual(values = mypal[3:18]) +
      facet_grid(~ group, scales = "free_x", switch = "x") +
      theme(
        axis.ticks.x = element_blank(),
        legend.position = "top",
        axis.text.x = element_blank(),
        strip.background = element_blank()
      ) +
      xlab("Groups") +
      scale_y_continuous(
        name = "Relative abundance (%)",
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25),
        labels = paste(seq(0, 100, 25))
      ) +
      labs(x = NULL, fill = name[1]) +
      theme(
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 12, colour = "black", family = "Arial"),
        axis.text = element_text(face = "bold", size = 12, color = "black", family = "Arial"),
        strip.text.x = element_text(face = "bold", size = 12, color = "black", family = "Arial"),
        panel.grid = element_blank(),
        legend.title = element_text(face = "bold", size = 12, color = "black", family = "Arial"),
        legend.text = element_text(face = "bold.italic", size = 12, color = "black", family = "Arial")
      )
    
    return(temp_plot)
  } else {
    print("Your ASV table is not equal.")
  }
}