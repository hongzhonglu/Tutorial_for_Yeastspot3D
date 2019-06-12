#----------------note
# this main script is used to handle with the gene mutation only from SNP information
# in this process, the gene with SNP will be translated into protein, based on which
# the SNP could be classified into nsSNP and sSNP
# Only nsSNP is used to mapping onto protein 3D structure

strain_classification <- read.table("data/strain_PDETOH_classification.txt", header = TRUE, stringsAsFactors = FALSE)
strain_classification <- strain_classification[, c('strain_name','type')]
strain_type0 <-"PDETOH_high"
strain_select1 <- chooseStrain(type = "PDETOH_high")

#------------batch process--------------------------------------------------------
# input the gene information
pdb_Ex <- read_excel("data/pdb_Ex refine for final residue distance calculation_manual check.xlsx")
pdb_Ex <- filter(pdb_Ex, is.na(pdb_Ex$With_distance))
pdb_Ex$pdbid <- paste(pdb_Ex$template, pdb_Ex$chain_new, sep = "@")
pdb_Ex <- select(pdb_Ex, locus, pdbid, qstart2, qend2, sstart2, send2)
geneWithSNP <- getGeneNameWithSNP()
pdb_Ex <- pdb_Ex[which(pdb_Ex$locus %in% geneWithSNP ==TRUE),]

#creat new file to store the results
outfile0 <- paste('result/hotspot from pdb_ex for ', strain_type0, sep = "")
dir.create(outfile0)
print(outfile0)

# start the batch process
for (i in 1:1047) {
  print(i)
  ss0 <- pdb_Ex$locus[i]
  mutated_gene0 <- preprocessSNP(ss0, gene_feature = gene_feature0)
  mutated_gene1 <- mutated_gene0[which(mutated_gene0$strain %in% strain_select1$Standardized_name), ]
  pdbID <- pdb_Ex$pdbid[i]
  distance_dir <- paste("residue_distance/pdb_ex/", pdbID, ".txt", sep = "")
  # run the function
  hotSpotAnalysis(
    gene0 = ss0,
    pdbID = pdbID,
    SNPlist0 = mutated_gene1,
    gene_annotation0 = gene_feature0,
    pdb = distance_dir,
    sstart0 = pdb_Ex$sstart2[i], # coordinate of orginal protein residues sequence
    send0 = pdb_Ex$send2[i],     # coordinate of orginal protein residues sequence
    qstart0 = pdb_Ex$qstart2[i], # coordinate of protein residues sequence in pdb file
    qend0 = pdb_Ex$qend2[i],     # coordinate of protein residues sequence in pdb file
    result_dir = outfile0,
    strain_type = strain_type0,
    input_dir=TRUE
  )
}



