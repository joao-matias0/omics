  ###########################Script do tratamento dos nossos dados##############
  
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
  
   #Replace Geneid column with CommonName
  #filtered_count_table <- filtered_count_table %>%
  #  mutate(Geneid = CommonName) %>%
  #  select(-CommonName)
  
  # Write the updated table to a file
  # write.table(filtered_count_table, "updated_count_table_with_names.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  #Remove rows where GENEID is NA or empty
  final_count_table <- filtered_count_table[!is.na(filtered_count_table$CommonName) & filtered_count_table$CommonName != "", ]
  final_count_table <- final_count_table[, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,19)]
  write.table(final_count_table, "final_count_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  #Generates new suffix names for repeating genes
  final_count_table[[ncol(final_count_table)]] <- make.names(final_count_table[[ncol(final_count_table)]], unique = TRUE)
  
  #Definir os nomes únicos como row names e remover a última coluna
  rownames(final_count_table) <- final_count_table[[ncol(final_count_table)]]
  final_count_table[[ncol(final_count_table)]] <- NULL
  ###########################Funtional Enrichment Analysis DESeq2##############
  ##############Differential_Expression_Analysis############################################
 library(DESeq2)
  # Step 1: Load the count table
  # Replace with your file path for the count data
  count_table <- read.table("filtered_count_table.txt", header = TRUE, sep = "\t", row.names = 1)
  
  # Ensure the count table contains numeric data only
  count_table <- as.matrix(count_table)
  count_table[count_table < 0] <- 0  # Replace negative values with zeros
  
  # Step 2: Define conditions
  # Adjust based on your experiment: first 6 columns are control, next 6 are drug1, last 6 are drug2
  conditions <- factor(c(rep("control", 6), rep("drug1", 6), rep("drug2", 6)))
  #Drug1= Methylone; Drug2= MDMA
  
  # Create metadata (colData) for DESeq2
  colData <- data.frame(condition = conditions, row.names = colnames(count_table))
  
  # Step 3: Create a DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_table,
                                colData = colData,
                                design = ~ condition)
  
  # Step 4: Run differential expression analysis
  dds <- DESeq(dds)
  
  # Step 5: Get results for Drug1 vs Control
  res_drug1 <- results(dds, contrast = c("condition", "drug1", "control"))
  
  # Get results for Drug2 vs Control
  res_drug2 <- results(dds, contrast = c("condition", "drug2", "control"))
  
  # Get results for Drug1 vs Drug2
  res_drug1VSdrug2 <- results(dds, contrast = c("condition", "drug1", "drug2"))
  
  # Step 6: Filter significant genes (adjusted p-value < 0.05)
  sig_genes_drug1 <- subset(res_drug1, padj < 0.05)
  sig_genes_drug2 <- subset(res_drug2, padj < 0.05)
  sig_genes_drug1VSdrug2 <- subset(res_drug1VSdrug2, padj < 0.05)
  
  # Step 7: Separate Upregulated and Downregulated Genes
  # Methylone
  upregulated_drug1 <- subset(sig_genes_drug1, log2FoldChange > 0)
  downregulated_drug1 <- subset(sig_genes_drug1, log2FoldChange < 0)
  
  # MDMA
  upregulated_drug2 <- subset(sig_genes_drug2, log2FoldChange > 0)
  downregulated_drug2 <- subset(sig_genes_drug2, log2FoldChange < 0)
  
  # Methylone vs MDMA
  upregulated_drug1VSDrug2 <- subset(sig_genes_drug1VSdrug2, log2FoldChange > 0)
  downregulated_drug1VSDrug2 <- subset(sig_genes_drug1VSdrug2, log2FoldChange < 0)
  
  # Step 8: Save Results to Files
  write.table(as.data.frame(upregulated_drug1), "upregulated_genes_methylone.txt", sep = "\t", row.names = TRUE, quote = FALSE)
  write.table(as.data.frame(downregulated_drug1), "downregulated_genes_methylone.txt", sep = "\t", row.names = TRUE, quote = FALSE)
  
  write.table(as.data.frame(upregulated_drug2), "upregulated_genes_mdma.txt", sep = "\t", row.names = TRUE, quote = FALSE)
  write.table(as.data.frame(downregulated_drug2), "downregulated_genes_mdma.txt", sep = "\t", row.names = TRUE, quote = FALSE)
  
  write.table(as.data.frame(upregulated_drug1VSDrug2), "upregulated_genes_methyloneVSmdma.txt", sep = "\t", row.names = TRUE, quote = FALSE)
  write.table(as.data.frame(downregulated_drug1VSDrug2), "downregulated_genes_methyloneVSmdma.txt", sep = "\t", row.names = TRUE, quote = FALSE)
  
  # Step 9: Print Summary
  cat("Drug1: Upregulated Genes:", nrow(upregulated_drug1), "\n")
  cat("Drug1: Downregulated Genes:", nrow(downregulated_drug1), "\n")
  
  cat("Drug2: Upregulated Genes:", nrow(upregulated_drug2), "\n")
  cat("Drug2: Downregulated Genes:", nrow(downregulated_drug2), "\n")
  
  cat("Drug1 vs Drug2: Upregulated Genes:", nrow(upregulated_drug1VSDrug2), "\n")
  cat("Drug1 vs Drug2: Downregulated Genes:", nrow(downregulated_drug1VSDrug2), "\n")

  up_reg_methylone <- read.table("upregulated_genes_methylone.txt") 
  up_reg_mdma <- read.table("upregulated_genes_mdma.txt", header = TRUE, sep = "\t") 
  down_reg_methylone <- read.table("downregulated_genes_methylone.txt", header = TRUE, sep = "\t") 
  down_reg_mdma <- read.table("downregulated_genes_mdma.txt", header = TRUE, sep = "\t")
  
  ##############Pathway analysis###############################################
  # Load required packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  #BiocManager::install(c("clusterProfiler", "org.Rn.eg.db"))
  
  # Load libraries
  library(clusterProfiler)
  library(org.Rn.eg.db)
  
  # Step 1: Load Gene Lists
  # Define function to read gene lists and extract IDs
  define_gene_list <- function(file_path) {
    gene_data <- read.table(file_path, header = TRUE, sep = "\t", skip = 1) # Skip the header
    return(gene_data[[1]]) # Extract the first column containing gene IDs
  }
  
  # Read gene lists
  upregulated_mdma <- define_gene_list("upregulated_genes_mdma.txt")
  downregulated_mdma <- define_gene_list("downregulated_genes_mdma.txt")
  upregulated_methylone <- define_gene_list("upregulated_genes_methylone.txt")
  downregulated_methylone <- define_gene_list("downregulated_genes_methylone.txt")
  
  # Step 2: Convert Ensembl IDs to Entrez IDs
  convert_ids <- function(ensembl_ids) {
    bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
  }
  
  upregulated_mdma_entrez <- convert_ids(upregulated_mdma)$ENTREZID
  downregulated_mdma_entrez <- convert_ids(downregulated_mdma)$ENTREZID
  upregulated_methylone_entrez <- convert_ids(upregulated_methylone)$ENTREZID
  downregulated_methylone_entrez <- convert_ids(downregulated_methylone)$ENTREZID
  
  # Step 3: Perform Enrichment Analysis
  perform_kegg <- function(entrez_ids) {
    enrichKEGG(
      gene = entrez_ids,
      organism = "rno",
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }
  
  perform_go <- function(entrez_ids) {
    enrichGO(
      gene = entrez_ids,
      OrgDb = org.Rn.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }
  
  # Run KEGG enrichment
  kegg_up_mdma <- perform_kegg(upregulated_mdma_entrez)
  kegg_down_mdma <- perform_kegg(downregulated_mdma_entrez)
  kegg_up_methylone <- perform_kegg(upregulated_methylone_entrez)
  kegg_down_methylone <- perform_kegg(downregulated_methylone_entrez)
  
  # Run GO enrichment
  go_up_mdma <- perform_go(upregulated_mdma_entrez)
  go_down_mdma <- perform_go(downregulated_mdma_entrez)
  go_up_methylone <- perform_go(upregulated_methylone_entrez)
  go_down_methylone <- perform_go(downregulated_methylone_entrez)
  
  # Step 4: Save Results
  save_results <- function(enrichment_result, file_name) {
    if (!is.null(enrichment_result)) {
      write.csv(as.data.frame(enrichment_result)[, c("ID", "Description", "GeneRatio", "p.adjust")], 
                file = file_name, row.names = FALSE)
    } else {
      cat("No significant pathways found for", file_name, "\n")
    }
  }
  
  # Save KEGG results | Bad results!!
  save_results(kegg_up_mdma, "kegg_up_mdma_results.csv")
  save_results(kegg_down_mdma, "kegg_down_mdma_results.csv")
  save_results(kegg_up_methylone, "kegg_up_methylone_results.csv")
  save_results(kegg_down_methylone, "kegg_down_methylone_results.csv")
  
  # Save GO results
  save_results(go_up_mdma, "go_up_mdma_results.csv")
  save_results(go_down_mdma, "go_down_mdma_results.csv")
  save_results(go_up_methylone, "go_up_methylone_results.csv")
  save_results(go_down_methylone, "go_down_methylone_results.csv")
  
  ##############Volcano Plots##################################################
  #Identifying Pathways
  #BiocManager::install("clusterProfiler")
  #BiocManager::install("enrichplot")
  
  #load libraries
  library(clusterProfiler)
  library(org.Rn.eg.db) 
  library(enrichplot)
  library(dplyr)
  
  # Install required packages if not already installed
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  # Uncomment these lines to install missing packages if needed
  # BiocManager::install("DESeq2")
  # BiocManager::install("ggplot2")
  # BiocManager::install("ggrepel")
  # BiocManager::install("org.Rn.eg.db")
  
  # Load libraries
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  
  # Step 1: Load the count table
  # Replace with your file path for the count data
  count_table <- read.table("filtered_count_table.txt", header = TRUE, sep = "\t", row.names = 1)
  
  
  # Ensure the count table contains numeric data only
  count_table <- as.matrix(count_table)
  count_table[count_table < 0] <- 0  # Replace negative values with zeros
  
  # Step 2: Define experimental conditions
  # Adjust based on your experiment: first 6 columns are control, next 6 are drug1, last 6 are drug2
  conditions <- factor(c(rep("control", 6), rep("drug1", 6), rep("drug2", 6)))
  
  # Create metadata (colData) for DESeq2
  colData <- data.frame(condition = conditions, row.names = colnames(count_table))
  
  # Step 3: Create a DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_table,
                                colData = colData,
                                design = ~ condition)
  
  # Step 4: Run differential expression analysis
  dds <- DESeq(dds)
  
  # Function to create volcano plots
  create_volcano_plot <- function(results_df, comparison, output_file) {
    # Convert DESeqResults to data frame
    results_df <- as.data.frame(results_df)
    
    # Add gene IDs
    results_df$Gene <- rownames(results_df)
    
    # Remove rows with NA values in log2FoldChange or padj
    results_df <- results_df[!is.na(results_df$log2FoldChange) & !is.na(results_df$padj), ]
    
    # Categorize genes for coloring
    results_df$Category <- "Not significant"  # Default category
    results_df$Category[results_df$padj < 0.05 & results_df$log2FoldChange > 0] <- "Upregulated"
    results_df$Category[results_df$padj < 0.05 & results_df$log2FoldChange < 0] <- "Downregulated"
    
    # Filter significant genes for labeling (top 30 by |log2FoldChange|)
    top_genes <- results_df %>%
      filter(Category %in% c("Upregulated", "Downregulated")) %>%
      arrange(-abs(log2FoldChange)) %>%
      head(30)
    
    # Map Ensembl IDs to common gene names (symbols)
    common_names <- AnnotationDbi::mapIds(org.Rn.eg.db,
                                          keys = top_genes$Gene,
                                          column = "SYMBOL",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")
    
    # Add common names to the top_genes data frame
    top_genes$GeneName <- common_names
    
    # If no common name exists, fall back to Ensembl ID
    top_genes$GeneName[is.na(top_genes$GeneName)] <- top_genes$Gene[is.na(top_genes$GeneName)]
    
    # Create the volcano plot
    volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = Category)) +
      geom_point(alpha = 0.8, size = 1.5) +  # Add points
      geom_text_repel(data = top_genes, 
                      aes(label = GeneName), 
                      size = 3, 
                      box.padding = 0.5, 
                      point.padding = 0.5) +  # Add labels for top genes
      scale_color_manual(values = c("Upregulated" = "red", 
                                    "Downregulated" = "blue", 
                                    "Not significant" = "gray")) +
      theme_minimal() +
      labs(title = paste("Volcano Plot:", comparison),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-Value") +
      theme(legend.title = element_blank())
    
    # Save the volcano plot to a file
    ggsave(output_file, plot = volcano_plot, width = 8, height = 6, dpi = 300)
    
    # Display the plot
    print(volcano_plot)
  }
  
  # Step 5: Define comparisons
  comparisons <- list(
    list(name = "Methylone vs Control", contrast = c("condition", "drug1", "control"), output = "volcano_plot_drug1_vs_control.png"),
    list(name = "MDMA vs Control", contrast = c("condition", "drug2", "control"), output = "volcano_plot_drug2_vs_control.png"),
    list(name = "Methylone vs MDMA", contrast = c("condition", "drug1", "drug2"), output = "volcano_plot_drug1_vs_drug2.png")
  )
  
  # Step 6: Run all comparisons and generate volcano plots
  for (comparison in comparisons) {
    results_df <- results(dds, contrast = comparison$contrast)
    create_volcano_plot(results_df, comparison$name, comparison$output)
  }
  
  
  
  ###########################Analysis
  ##############Common_Pathways_GO#############################################
  ##DownReg
  # Load necessary libraries
  library(dplyr)
  
  # Load the GO term data
  go_mdma <- read.csv("go_down_mdma_results.csv")
  go_methylone <- read.csv("go_down_methylone_results.csv")
  
  # Ensure the GO term column names are consistent
  # Replace "GO.Term" with the actual column name if it's different
  colnames(go_mdma) <- make.names(colnames(go_mdma))
  colnames(go_methylone) <- make.names(colnames(go_methylone))
  
  # Find the intersection of GO terms
  common_go_terms <- intersect(go_mdma$ID, go_methylone$GO.Term)
  
  # Save the common GO terms to a CSV file
  write.csv(common_go_terms, "common_go_terms_down.csv", row.names = FALSE)
  
  # Print the results
  print("Common GO terms found:")
  print(common_go_terms)
  
  ##UpReg
  # Load the GO term data
  go_mdma <- read.csv("go_up_mdma_results.csv")
  go_methylone <- read.csv("go_up_methylone_results.csv")
  
  # Ensure the GO term column names are consistent
  # Replace "GO.Term" with the actual column name if it's different
  colnames(go_mdma) <- make.names(colnames(go_mdma))
  colnames(go_methylone) <- make.names(colnames(go_methylone))
  
  # Find the intersection of GO terms
  common_go_terms <- intersect(go_mdma$ID, go_methylone$GO.Term)
  
  # Save the common GO terms to a CSV file
  write.csv(common_go_terms, "common_go_terms_up.csv", row.names = FALSE)
  
  # Print the results
  print("Common GO terms found:")
  print(common_go_terms)
  
  ##############Common_Genes###################################################
  ##DownReg
  # Load necessary libraries
  library(dplyr)
  library(tibble) #needed for rownames_to_collumn
  library(org.Rn.eg.db)
  
  # Load the datasets
  down_gene_mdma <- read.table("downregulated_genes_mdma.txt", header = TRUE, row.names = 1)
  down_gene_methylone <- read.table("downregulated_genes_methylone.txt", header = TRUE, row.names = 1)
  
  # Convert row names to a column named "GeneID"
  down_gene_mdma <- down_gene_mdma %>%
    rownames_to_column(var = "GeneID")
  
  down_gene_methylone <- down_gene_methylone %>%
    rownames_to_column(var = "GeneID")
  
  # Find the common GeneIDs
  common_genes_down <- intersect(down_gene_mdma$GeneID, down_gene_methylone$GeneID)
  
  # Save the common genes to a CSV file
  write.csv(common_genes_down, "common_genes_down.csv", row.names = FALSE)
  
  # Print the results
  print("Common GeneIDs found:")
  print(common_genes_down)
  
  ##UpReg
  # Load the datasets
  up_gene_mdma <- read.table("upregulated_genes_mdma.txt", header = TRUE, row.names = 1)
  up_gene_methylone <- read.table("upregulated_genes_methylone.txt", header = TRUE, row.names = 1)
  
  # Convert row names to a column named "GeneID"
  up_gene_mdma <- up_gene_mdma %>%
    rownames_to_column(var = "GeneID")
  
  up_gene_methylone <- up_gene_methylone %>%
    rownames_to_column(var = "GeneID")
  
  # Find the common GeneIDs
  common_genes_up <- intersect(up_gene_mdma$GeneID, up_gene_methylone$GeneID)
  
  # Save the common genes to a CSV file
  write.csv(common_genes_up, "common_genes_up.csv", row.names = FALSE)
  
  # Print the results
  print("Common GeneIDs found:")
  print(common_genes_up)
  
  ####
  ####
  ##Transform into common names
  
  common_names_gene_down <- read.table("common_genes_down.csv", header = TRUE)
  common_names_gene_up <- read.table("common_genes_up.csv", header = TRUE)
  
  common_names_down <- AnnotationDbi::mapIds(
    org.Rn.eg.db,
    keys = common_names_gene_down$x,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  ) 
  
  common_names_up <- AnnotationDbi::mapIds(
    org.Rn.eg.db,
    keys = common_names_gene_up$x,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  ) 
  
  write.csv(common_names_down, "common_genes_down_names.csv", row.names = FALSE)
  write.csv(common_names_up, "common_genes_up_names.csv", row.names = FALSE)
  
  
  
  ##############Heatmaps#######################################################
  ##############Top10_Enriched_terms#######################################################
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  
  ##Upregulated - MDMA
  # Load the dataset
  go_data_up_mdma <- read.csv("go_up_mdma_results.csv", header = TRUE)
  
  # Ensure the columns are properly named
  # Replace "GO.Term", "Description", and "P.value" with actual column names if different
  colnames(go_data_up_mdma) <- make.names(colnames(go_data_up_mdma))
  
  # Create a new column for -log10(P-value)
  go_data_up_mdma <- go_data_up_mdma %>%
    mutate(`-log10(P)` = -log10(p.adjust)) %>%
    arrange(desc(`-log10(P)`))  # Sort by significance
  
  # Select the top 10 enriched terms
  top10_go_up_mdma <- head(go_data_up_mdma, 10)
  write.csv(top10_go_up_mdma, "top10_enriched_go_up_terms_mdma.csv", row.names = FALSE)
  
  #Upregulated - Methylone
  
  go_data_up_methylone <- read.csv("go_up_methylone_results.csv", header = TRUE)
  
  # Ensure the columns are properly named
  # Replace "GO.Term", "Description", and "P.value" with actual column names if different
  colnames(go_data_up_methylone) <- make.names(colnames(go_data_up_methylone))
  
  # Create a new column for -log10(P-value)
  go_data_up_methylone <- go_data_up_methylone %>%
    mutate(`-log10(P)` = -log10(p.adjust)) %>%
    arrange(desc(`-log10(P)`))  # Sort by significance
  
  # Select the top 10 enriched terms
  top10_go_up_methylone <- head(go_data_up_methylone, 10)
  write.csv(top10_go_up_methylone, "top10_enriched_go_up_terms_methylone.csv", row.names = FALSE)

  ##Downregulated - MDMA
  # Load the dataset
  go_data_down_mdma <- read.csv("go_down_mdma_results.csv", header = TRUE)
  
  # Ensure the columns are properly named
  # Replace "GO.Term", "Description", and "P.value" with actual column names if different
  colnames(go_data_down_mdma) <- make.names(colnames(go_data_down_mdma))
  
  # Create a new column for -log10(P-value)
  go_data_down_mdma <- go_data_down_mdma %>%
    mutate(`-log10(P)` = -log10(p.adjust)) %>%
    arrange(desc(`-log10(P)`))  # Sort by significance
  
  # Select the top 10 enriched terms
  top10_go_down_mdma <- head(go_data_down_mdma, 10)
  write.csv(top10_go_down_mdma, "top10_enriched_go_down_terms_mdma.csv", row.names = FALSE)
  
  
  #Downregulated - Methylone
  # Load the dataset
  go_data_down_methylone <- read.csv("go_down_methylone_results.csv", header = TRUE)
  
  # Ensure the columns are properly named
  # Replace "GO.Term", "Description", and "P.value" with actual column names if different
  colnames(go_data_down_methylone) <- make.names(colnames(go_data_down_methylone))
  
  # Create a new column for -log10(P-value)
  go_data_down_methylone <- go_data_down_methylone %>%
    mutate(`-log10(P)` = -log10(p.adjust)) %>%
    arrange(desc(`-log10(P)`))  # Sort by significance
  
  # Select the top 10 enriched terms
  top10_go_down_methylone <- head(go_data_down_methylone, 10)
  write.csv(top10_go_down_methylone, "top10_enriched_go_down_terms_methylone.csv", row.names = FALSE)
  
  
  
  
  ##############Volcano_Plots_2_TEST###########################################
  
  # Identifying Pathways
  # Uncomment and install packages if necessary
  # BiocManager::install("clusterProfiler")
  # BiocManager::install("enrichplot")
  
  # Load libraries
  library(clusterProfiler)
  library(org.Rn.eg.db)
  library(enrichplot)
  library(dplyr)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationDbi)
  
  # Install required packages if not already installed
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  # Load count table
  # Replace with your file path for the count data
  count_table <- read.table("filtered_count_table.txt", header = TRUE, sep = "\t", row.names = 1)
  
  # Ensure the count table contains numeric data only
  count_table <- as.matrix(count_table)
  count_table[count_table < 0] <- 0  # Replace negative values with zeros
  
  # Define experimental conditions
  # Adjust based on your experiment: first 6 columns are control, next 6 are drug1, last 6 are drug2
  conditions <- factor(c(rep("control", 6), rep("drug1", 6), rep("drug2", 6)))
  
  # Create metadata (colData) for DESeq2
  colData <- data.frame(condition = conditions, row.names = colnames(count_table))
  
  # Create a DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_table,
                                colData = colData,
                                design = ~ condition)
  
  # Run differential expression analysis
  dds <- DESeq(dds)
  
  # Function to create volcano plots
  create_volcano_plot <- function(results_df, comparison, output_file, highlight_genes = NULL) {
    # Convert DESeqResults to data frame
    results_df <- as.data.frame(results_df)
    
    # Add gene IDs
    results_df$Gene <- rownames(results_df)
    
    # Remove rows with NA values in log2FoldChange or padj
    results_df <- results_df[!is.na(results_df$log2FoldChange) & !is.na(results_df$padj), ]
    
    # Categorize genes for coloring
    results_df$Category <- "Not significant"  # Default category
    results_df$Category[results_df$padj < 0.05 & results_df$log2FoldChange > 0] <- "Upregulated"
    results_df$Category[results_df$padj < 0.05 & results_df$log2FoldChange < 0] <- "Downregulated"
    
    # Highlight genes if provided
    if (!is.null(highlight_genes)) {
      # Map common names to GeneID format
      highlight_gene_ids <- AnnotationDbi::mapIds(org.Rn.eg.db,
                                                  keys = highlight_genes,
                                                  column = "ENSEMBL",
                                                  keytype = "SYMBOL",
                                                  multiVals = "first")
      
      # Ensure matching GeneIDs are identified
      results_df$Highlight <- ifelse(results_df$Gene %in% highlight_gene_ids, "Highlight", "No Highlight")
    } else {
      results_df$Highlight <- "No Highlight"
    }
    
    # Filter significant genes for labeling (top 30 by |log2FoldChange|)
    top_genes <- results_df %>%
      filter(Category %in% c("Upregulated", "Downregulated")) %>%
      arrange(-abs(log2FoldChange)) %>%
      head(30)
    
    # Map Ensembl IDs to common gene names (symbols)
    common_names <- AnnotationDbi::mapIds(org.Rn.eg.db,
                                          keys = top_genes$Gene,
                                          column = "SYMBOL",
                                          keytype = "ENSEMBL",
                                          multiVals = "first")
    
    # Add common names to the top_genes data frame
    top_genes$GeneName <- common_names
    
    # If no common name exists, fall back to Ensembl ID
    top_genes$GeneName[is.na(top_genes$GeneName)] <- top_genes$Gene[is.na(top_genes$GeneName)]
    
    # Create the volcano plot
    volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = Category)) +
      geom_point(alpha = 0.8, size = 1.5) +  # Add points
      geom_point(data = subset(results_df, Highlight == "Highlight"), 
                 aes(x = log2FoldChange, y = -log10(padj)), 
                 color = "green", size = 2, shape = 21, stroke = 1.5) + # Highlight specific genes
      geom_text_repel(data = top_genes, 
                      aes(label = GeneName), 
                      size = 3, 
                      box.padding = 0.5, 
                      point.padding = 0.5) +  # Add labels for top genes
      scale_color_manual(values = c("Upregulated" = "red", 
                                    "Downregulated" = "blue", 
                                    "Not significant" = "gray")) +
      theme_minimal() +
      labs(title = paste("Volcano Plot:", comparison),
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-Value") +
      theme(legend.title = element_blank())
    
    # Save the volcano plot to a file
    ggsave(output_file, plot = volcano_plot, width = 8, height = 6, dpi = 300)
    
    # Display the plot
    print(volcano_plot)
  }
  
  # Read the gene list file
  highlight_genes <- readLines("highlight_genes.txt")
  
  # Define comparisons
  comparisons <- list(
    list(name = "Methylone vs Control", contrast = c("condition", "drug1", "control"), output = "volcano_plot_drug1_vs_control.png"),
    list(name = "MDMA vs Control", contrast = c("condition", "drug2", "control"), output = "volcano_plot_drug2_vs_control.png"),
    list(name = "Methylone vs MDMA", contrast = c("condition", "drug1", "drug2"), output = "volcano_plot_drug1_vs_drug2.png")
  )
  
  # Run all comparisons and generate volcano plots
  for (comparison in comparisons) {
    results_df <- results(dds, contrast = comparison$contrast)
    create_volcano_plot(results_df, comparison$name, comparison$output, highlight_genes)
  }
  