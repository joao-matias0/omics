#install necessary packages (only need to run once)
install.packages("BiocManager")
install.packages("readr")
BiocManager::install('clusterProfiler')
BiocManager::install('DOSE')
BiocManager::install('org.Hs.eg.db')

#Transform a .txt table into a .csv
data <- read.table("final_count_table.txt", header = TRUE, sep = "\t")
write.csv(data, "final_count_table.txt", row.names = FALSE)


























#load table with list of cancer driver neighbours 

library(readr)
neib <- read_csv("neibrho.csv")

signeibplus=neib$Group.1[neib$z>=3.1 & neib$sign=="positive" & neib$psurvpadj<=0.05]
univgenes=neib$Group.1

#The $ sign extracts information for the given collumn | [neib$z>=3.1] extracts
#the rows that have a z value equal or above 3.1

#perform over-representation analysis of GO BP terms

library(clusterProfiler)
library(org.Hs.eg.db)

ego <- enrichGO(gene          = signeibplus,
                universe      = univgenes,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "BP",   #BP:: Biological Process (Gene Ontology)
                pAdjustMethod = "BH",   #Method for multiple correction BH Benjamin Haufer controls false discovery rate
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 1,
                minGSSize = 20,         # 
                maxGSSize = 250,        #
                readable      = TRUE)

library(DOSE)

#    ego@result$FoldEnrichment=parse_ratio(ego@result$GeneRatio)/parse_ratio(ego@result$BgRatio)
#Not needed anymore as the newer package version already does this automatically

#to check the significant results in a table
egodf=as.data.frame(ego)

#to remove enriched terms with small number of associated genes
todrop=ego@result$ID[ego@result$Count<20]
ego=dropGO(ego,term=todrop)

egodf=as.data.frame(ego)

#to remove redundancy

egosimp=simplify(
    ego,
    cutoff = 0.8,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = NULL
)

egosimpdf=as.data.frame(egosimp)

#to view results in a plot

dotplot(egosimp, x="FoldEnrichment",color="p.adjust",size="Count",showCategory=20,orderBy="Count")

#to export the results table

write.table(egodf,"egodf.txt",row.names=F)

#perform over-representation analysis of KEGG pathways

#this analysis requires using a different gene id
gene.df <- bitr(gene= signeibplus, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

univ.df <- bitr(gene= univgenes, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene          = gene.df$ENTREZID,
                 organism     = 'hsa',                 #Homo Sapiens
                 pAdjustMethod = "BH",
                 universe= univ.df$ENTREZID,
                 pvalueCutoff = 0.05)

kk@result$FoldEnrichment=parse_ratio(kk@result$GeneRatio)/parse_ratio(kk@result$BgRatio)

kkdf=as.data.frame(kk)

dotplot(kk, x="FoldEnrichment",color="p.adjust",size="Count",showCategory=20,orderBy="Count")

#view results on kegg map
browseKEGG(kk, 'hsa03460')


