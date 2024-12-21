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
heatmap(heatmap_data, 
        main = "Heatmap of First 6 Columns", 
        xlab = "Samples", 
        ylab = "Genes", 
        col = heat.colors(256), 
        scale = "row")


