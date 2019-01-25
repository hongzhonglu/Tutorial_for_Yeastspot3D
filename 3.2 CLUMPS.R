#----------------note
#this main script is used to handle with the gene mutation only from SNP information
#in this process, the gene with SNP will be translated into protein, based on which
#the SNP could be classified into nsSNP and sSNP
#Only nsSNP is used to mapping onto protein 3D structure

#source('genome_annotation.R')

#step0 choose samples that need to be analyzed
strain_classification <- read_excel("data/strain_classification.xls")
unique(strain_classification$Clades)

strain_type <- "bioethanol"
strain_select1 <- filter(strain_classification, str_detect(strain_classification$Clades, strain_type)) %>%
  select(., Standardized_name)

##batch process
# input the gene information
pdb_Ex <- read_excel("data/pdb_Ex refine for final residue distance calculation_manual check.xlsx")
pdb_Ex <- filter(pdb_Ex, is.na(pdb_Ex$With_distance))
pdb_Ex$pdbid <- paste(pdb_Ex$template, pdb_Ex$chain_new, sep = "@")
pdb_Ex <- select(pdb_Ex, locus, pdbid, qstart2, qend2, sstart2, send2)
geneWithSNP <- getGeneNameWithSNP()
pdb_Ex <- pdb_Ex[which(pdb_Ex$locus %in% geneWithSNP ==TRUE),]

#add two more clumns
pdb_Ex$strain_type <- strain_type
pdb_Ex$p_value <- NA

#creat new file to store the results
outfile0 <- paste('result/CLUMPS from pdb_ex for ', strain_type, sep = "")
dir.create(outfile0)
print(outfile0)


#main loop
for (i in 1:1047) {
  print(i)
  i <- 1
  ss0 <- pdb_Ex$locus[i]
  mutated_gene0 <- preprocessSNP(ss0)
  mutated_gene1 <- mutated_gene0[which(mutated_gene0$strain %in% strain_select1$Standardized_name), ]
  pdbID <- pdb_Ex$pdbid[i]
  p1 <- pdb_Ex$sstart2[i]
  p2 <- pdb_Ex$send2[i]
  #in some function, gene_snp was not put as the input, should be careful
  #gene_snp <- getGeneCoordinate(gene_name = ss0, genesum = gene_feature_GEM)
  result0 <- clumpsAnalysis(gene0 = ss0, SNPlist0 = mutated_gene1, gene_annotation0 = gene_feature_GEM, pdbID0 = pdbID, sstart0 = p1, send0 = p2)
  pdb_Ex$p_value[i] <- result0
}

# save the result
write.table(pdb_Ex, paste(outfile0,'/','pdb_EX.txt', sep = ""), row.names = FALSE, sep = "\t")


