library(gplots)

library(org.Ce.eg.db)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(enrichplot)

#Use Wormbase IDs to create a character vector of gene set of interest
genes <- Brood_sig #set of 448 genes
genes2<-NoiseGenes #set of 114 genes
genes3<-ParentGenes #set of 97 genes
genes4<-UpGenes #positively associated 159 genes out of 448
genes5<-DownGenes #negatively associated 289 out of 448
genes6<-Beta1Beta2_ELO #universe background set of 8824 genes
genes7<-EnrichedGenes_4 #set of 246 genes from 10-gene prediction


# Suppose your gene list is stored in a character vector called 'genes'
# Convert your gene list to official gene symbols
gene_symbols <- mapIds(org.Ce.eg.db, keys = genes, keytype = "WORMBASE", column = "SYMBOL", multiVals = "first")
gene_symbols2 <- mapIds(org.Ce.eg.db, keys = genes2, keytype = "WORMBASE", column = "SYMBOL", multiVals = "first")
gene_symbols3 <- mapIds(org.Ce.eg.db, keys = genes3, keytype = "WORMBASE", column = "SYMBOL", multiVals = "first")
gene_symbols4 <- mapIds(org.Ce.eg.db, keys = genes4, keytype = "WORMBASE", column = "SYMBOL", multiVals = "first")
gene_symbols5 <- mapIds(org.Ce.eg.db, keys = genes5, keytype = "WORMBASE", column = "SYMBOL", multiVals = "first")
gene_symbols6 <- mapIds(org.Ce.eg.db, keys = genes8, keytype = "WORMBASE", column = "SYMBOL", multiVals = "first")
gene_symbols7 <- mapIds(org.Ce.eg.db, keys = genes9, keytype = "WORMBASE", column = "SYMBOL", multiVals = "first")


ego <- enrichGO(gene = gene_symbols, OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego2 <- enrichGO(gene = gene_symbols2, OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego3 <- enrichGO(gene = gene_symbols3, OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego4 <- enrichGO(gene = gene_symbols4, OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego5 <- enrichGO(gene = gene_symbols5, OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego6 <- enrichGO(gene = gene_symbols, universe = gene_symbols6,OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego7 <- enrichGO(gene = gene_symbols2, universe = gene_symbols6,OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego8 <- enrichGO(gene = gene_symbols3, universe = gene_symbols6,OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego9 <- enrichGO(gene = gene_symbols7, OrgDb = org.Ce.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
#ego6,7,8 above include the 8824-gene background set as the universe


?enrichGO
#example of how to write table given ego object
write.table(ego@result,"GOterms_448genes_check.txt",quote = F,sep = "\t")


# Plot the enriched GO terms
dotplot(ego, showCategory=20)


