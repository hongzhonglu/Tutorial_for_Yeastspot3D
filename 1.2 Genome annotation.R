#gene feature summary
gene_feature0 <- mergeGeneAnnotationFeature()
gene_feature0 <- qualityCheckFromCDStoProtein(gene_feature0)
write.table(gene_feature0,'result/gene_feature0.txt', row.names = FALSE, sep = "\t")
