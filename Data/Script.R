
  #Script do tratamento dos nossos dados
  
  file1 <- read.table("1-FC_counts.txt", header = TRUE, sep = "\t")
  file2 <- read.table("2-FC_counts.txt", header = TRUE, sep = "\t")
  file3 <- read.table("3-FC_counts.txt", header = TRUE, sep = "\t")
  file4 <- read.table("4-FC_counts.txt", header = TRUE, sep = "\t")
  file5 <- read.table("5-FC_counts.txt", header = TRUE, sep = "\t")
  file6 <- read.table("6-FC_counts.txt", header = TRUE, sep = "\t")
  file7 <- read.table("7-FC_counts.txt", header = TRUE, sep = "\t")
  file8 <- read.table("8-FC_counts.txt", header = TRUE, sep = "\t")
  file9 <- read.table("9-FC_counts.txt", header = TRUE, sep = "\t")
  file10 <- read.table("10-FC_counts.txt", header = TRUE, sep = "\t")
  file11 <- read.table("11-FC_counts.txt", header = TRUE, sep = "\t")
  file12 <- read.table("12-FC_counts.txt", header = TRUE, sep = "\t")
  file13 <- read.table("13-FC_counts.txt", header = TRUE, sep = "\t")
  file14 <- read.table("14-FC_counts.txt", header = TRUE, sep = "\t")
  file15 <- read.table("15-FC_counts.txt", header = TRUE, sep = "\t")
  file16 <- read.table("16-FC_counts.txt", header = TRUE, sep = "\t")
  file17 <- read.table("17-FC_counts.txt", header = TRUE, sep = "\t")
  file18 <- read.table("18-FC_counts.txt", header = TRUE, sep = "\t")
  
  # Create a list of all the files
  file_list <- list(file1, file2, file3, file4, file5, file6, file7, file8, file9, 
                    file10, file11, file12, file13, file14, file15, file16, file17, file18)
  
  # Use Reduce to merge all the files by "Geneid"
  count_table <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), file_list)
  
  #Create a file
  write.table(count_table, "count_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  ###Removing rows where sum is 0
  # Load the table
  count_table <- read.table("count_table.txt", header = TRUE, sep = "\t")
  
  # Calculate row sums excluding the first column (Geneid)
  row_sums <- rowSums(count_table[ , -1])
  
  # Filter out rows where the sum is 0
  #filtered_table <- count_table[row_sums != 0, ]
  
  #Filter out rows where the sum is 10
  filtered_table <- count_table[row_sums >= 10, ]
  
  # Write the filtered table to a new file
  write.table(filtered_table, "filtered_count_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  
  #############Repalcing Geneid with Gene Names
  #Install BiocManager
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  #BiocManager::install("org.Rn.eg.db")
  
  # Load required libraries
  library(org.Rn.eg.db)
  library(dplyr)
  
  # Load the count table
  filtered_count_table <- read.table("filtered_count_table.txt", header = TRUE, sep = "\t")
  
  # Map Ensembl IDs to symbols (common names)
  # The mapIds function retrieves mappings from org.Rn.eg.db
  common_names <- AnnotationDbi::mapIds(
    org.Rn.eg.db,
    keys = count_table$Geneid,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Add common names as a new column
  filtered_count_table <- filtered_count_table %>%
    mutate(CommonName = common_names[Geneid])
  
  # Move the original Geneid column to the last position, if needed
  filtered_count_table <- filtered_count_table %>%
    relocate(Geneid, .after = last_col())
  
  # Replace Geneid column with CommonName
  #filtered_count_table <- filtered_count_table %>%
  #  mutate(Geneid = CommonName) %>%
  #  select(-CommonName)
  
  # Write the updated table to a file
  write.table(filtered_count_table, "updated_count_table_with_names.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  #Remove rows where GENEID is NA or empty
  final_count_table <- filtered_count_table[!is.na(filtered_count_table$CommonName) & filtered_count_table$CommonName != "", ]
  final_count_table <- final_count_table[, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,19)]
  write.table(final_count_table, "final_count_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  #Generates new suffix names for repeating genes
  final_count_table[[ncol(final_count_table)]] <- make.names(final_count_table[[ncol(final_count_table)]], unique = TRUE)
  
  #Definir os nomes únicos como row names e remover a última coluna
  rownames(final_count_table) <- final_count_table[[ncol(final_count_table)]]
  final_count_table[[ncol(final_count_table)]] <- NULL
  
  #Write final_count_table
  write.table(final_count_table, "Fire.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  #########################################################################################
  ####Funtional Enrichment Analysis
  
  #Defining object for the functonal enrichment
  data <- final_count_table
  
  #Identifying Pathways
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("enrichplot")
  
  #load libraries
  library(clusterProfiler)
  library(org.Rn.eg.db) 
  library(enrichplot)
  library(dplyr)
  
  
  
  # Subset the first 6 columns (excluding row names)
  heatmap_data <- as.matrix(data[, 1:6])
  
  # Optional: Normalize the data (e.g., scale rows)
  heatmap_data <- scale(heatmap_data)
  
  # Create a heatmap
  #heatmap(heatmap_data, 
         # main = "Heatmap of First 6 Columns", 
          #xlab = "Samples", 
         # ylab = "Genes", 
         # col = heat.colors(256), 
         # scale = "row")
  
  
  
  # Install required packages if not already installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
  if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2")
  if (!requireNamespace("ggrepel", quietly = TRUE)) BiocManager::install("ggrepel")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
  if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) BiocManager::install("org.Rn.eg.db")
  
  # Load libraries
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  
  # Step 1: Load the final count table
  final_count_table <- read.table("final_count_table.txt", header = TRUE, sep = "\t")
  
  # Step 2: Ensure numeric columns are selected for counts
  # Remove non-numeric columns (e.g., Geneid, CommonName) if present
  count_data <- final_count_table[, sapply(final_count_table, is.numeric)]
  
  # Convert to a matrix
  count_data <- as.matrix(count_data)
  
  # Ensure all values are numeric
  mode(count_data) <- "numeric"
  
  # Replace negative values with zeros (if any)
  count_data[count_data < 0] <- 0
  
  # Step 3: Define experimental conditions
  # Adjust based on your experiment: first 6 columns are control, next 6 are drug1, last 6 are drug2
  conditions <- factor(c(rep("control", 6), rep("drug1", 6), rep("drug2", 6)))
  
  # Create metadata (colData)
  colData <- data.frame(condition = conditions, row.names = colnames(count_data))
  
  # Step 4: Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = colData, design = ~ condition)
  
  # Step 5: Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Step 6: Define a function to create volcano plots
  create_volcano_plot <- function(results_df, comparison, output_file) {
    # Convert DESeqResults to a data frame
    results_df <- as.data.frame(results_df)
    results_df$Gene <- rownames(results_df)
    
    # Remove NA values in log2FoldChange and padj
    results_df <- results_df[!is.na(results_df$log2FoldChange) & !is.na(results_df$padj), ]
    
    # Categorize genes
    results_df$Category <- "Not significant"
    results_df$Category[results_df$padj < 0.05 & results_df$log2FoldChange > 0] <- "Upregulated"
    results_df$Category[results_df$padj < 0.05 & results_df$log2FoldChange < 0] <- "Downregulated"
    
    # Identify top genes for labeling
    top_genes <- results_df %>%
      filter(Category %in% c("Upregulated", "Downregulated")) %>%
      arrange(-abs(log2FoldChange)) %>%
      head(10)
    
    # Determine valid keytype for mapping
    valid_keytype <- "ENSEMBL"  # Default keytype; adjust based on your gene IDs
    if (!all(top_genes$Gene %in% keys(org.Rn.eg.db, keytype = valid_keytype))) {
      cat("Some gene keys are not valid for keytype:", valid_keytype, "\n")
      cat("Attempting to map using an alternative keytype: SYMBOL...\n")
      valid_keytype <- "SYMBOL"  # Switch to SYMBOL if ENSEMBL keys are not valid
    }
    
    # Map gene IDs to common names
    common_names <- AnnotationDbi::mapIds(org.Rn.eg.db,
                                          keys = top_genes$Gene,
                                          column = "SYMBOL",
                                          keytype = valid_keytype,
                                          multiVals = "first")
    
    # Add common names to top_genes
    top_genes$GeneName <- common_names
    top_genes$GeneName[is.na(top_genes$GeneName)] <- top_genes$Gene[is.na(top_genes$GeneName)]
    
    # Create the volcano plot
    volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = Category)) +
      geom_point(alpha = 0.8, size = 1.5) +
      geom_text_repel(data = top_genes, aes(label = GeneName), size = 3, box.padding = 0.5, point.padding = 0.5) +
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
      theme_minimal() +
      labs(title = paste("Volcano Plot:", comparison), x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
      theme(legend.title = element_blank())
    
    # Save the plot
    ggsave(output_file, plot = volcano_plot, width = 8, height = 6, dpi = 300)
    
    # Display the plot
    print(volcano_plot)
  }
  
  # Step 7: Define comparisons
  comparisons <- list(
    list(name = "Drug1 vs Control", contrast = c("condition", "drug1", "control"), output = "volcano_plot_drug1_vs_control.png"),
    list(name = "Drug2 vs Control", contrast = c("condition", "drug2", "control"), output = "volcano_plot_drug2_vs_control.png"),
    list(name = "Drug1 vs Drug2", contrast = c("condition", "drug1", "drug2"), output = "volcano_plot_drug1_vs_drug2.png")
  )
  
  # Step 8: Run all comparisons and generate volcano plots
  for (comparison in comparisons) {
    results_df <- results(dds, contrast = comparison$contrast)
    create_volcano_plot(results_df, comparison$name, comparison$output)
  }
  
  
  