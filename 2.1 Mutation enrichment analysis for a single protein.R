# this script is showing how to analyze the data for the SNP data in the excel format
# it is also base for hotspot and clumps analysis method
# the SNP format is as followed:
# Chr       pos      Ref       Alt
# chr1      19971    C         T
# chr1      19972    T         C


# function need to adjusted:
PositionResidueSNP <- function(mutatedPosition, alted, geneName, gene_feature) {
  #mutatedPosition = 130975
  #alted ='A'
  #geneName = "YAL012W"
  gene_snp <- getGeneCoordinate(gene_name = geneName, genesum = gene_feature)
  mutation_position <- which(gene_snp[['gene_coordinate']]==mutatedPosition)

  gene_snp[['gene']][mutation_position] <- alted

  # translation
  #library(seqinr)
  realcds <- str_to_lower(paste(gene_snp[["gene"]], collapse = ""))
  toycds <- s2c(realcds)
  gene_snp[["protein_mutated"]] <- translate(seq = toycds)

  # find the relative postion of mutated amino acids
  aa_position <- which(gene_snp[["protein"]] != gene_snp[["protein_mutated"]])
  aa_type <- gene_snp[["protein_mutated"]][aa_position]

  # built the relation between aa_position and aa_type
  # aa_postion and aa_type should contain one element
  mutatedAA <- paste(aa_type, aa_position, sep = "@@") # this estabolish the relation between the postion and mutated amino acids
  return(mutatedAA)
}


printSNPforGene <- function (gene0 = ss0, SNPlist0 = mutated_gene1, gene_annotation0 = gene_feature0,
                             pdbID0 = pdbID, sstart0 = p1, send0 = p2)
{
  ss <- gene0
  p3 <- paste(sstart0, send0, sep = "-")
  seq_3D_origin <- sstart0:send0
  pos_residue1 <- list()
  for (j in 1:nrow(SNPlist0)) {
    cat("Process all the SNPs from SNP list to obtain the SNP belong to the input gene:")
    print(j)
    pos_residue1[[j]] <- PositionResidueSNP(SNPlist0$Pos[j],
                                            SNPlist0$Alt[j], ss, gene_feature = gene_annotation0)
  }
  pos_residue_df <- ResidueSum(pos_residue1)
  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_annotation0)
  gene_snp[["pro_coordinate"]] <- 1:length(gene_snp[["protein"]])
  gene_snp[["residue"]] <- getMultipleReactionFormula(pos_residue_df$residue,
                                                      pos_residue_df$pos, gene_snp[["pro_coordinate"]])
  residue_3D <- gene_snp[["residue"]][seq_3D_origin]
  residue_3D0 <- residue_3D[!is.na(residue_3D)]
  residue_3D0 <- str_split(residue_3D0, ";")
  residue_3D0 <- unlist(residue_3D0)
  unique(residue_3D0)
  tmp <- table(residue_3D0)
  result <- as.data.frame(tmp)
  result <- result %>% separate(., residue_3D0, into = c("alt",
                                                         "position"), sep = "@@")
  result$position <- as.numeric(result$position)
  result <- result %>% arrange(., position)
  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_annotation0)
  result$ref <- getSingleReactionFormula(gene_snp[["protein"]],
                                         gene_snp[["protein_coordinate"]], result$position)
  result$orf <- ss
  result$pdbID <- pdbID0
  result <- select(result, orf, ref, position, alt, Freq, pdbID)
  cat("The SNP summary:")
  print(result)
  result0 <- select(result, orf, ref, position, alt)
  result0 <- result0 %>% unite(mutation, c("ref", "position",
                                           "alt"), sep = "")
  result$mutation <- result0$mutation
  return(result)
}


#get the gene name
#try to calculate the mutation on the amino acids based on the coordination on the chromosome
mutated_test <- read_excel("data/snp_adaption_to_high_ethanol.XLS")
mutated_test$Chr <- str_trim(mutated_test$Chr, side = "both")
mutated_test$Pos <- as.numeric(mutated_test$Pos)



for (i in seq(length(mutated_test$Chr))){
  print(i)
  mutated_test$Gene2[i] <- getGeneName(mutated_test$Chr[i],mutated_test$Pos[i])
}

mutated_test0 <- filter(mutated_test, Gene2 != "INTERGENIC") ##filter the mutated test

#choose the metabolic gene
#if the gene is type of "complement", then the complement_sign is "TRUE"
#else the complement_sign is "FALSE"
gene_feature0$complement_sign <- str_detect(gene_feature0$cds_location,"complement")
index_m <- which(mutated_test0$Gene2 %in% gene_feature0$locus_tag ==TRUE)
mutated_gene <- mutated_test0[index_m,]


mutated_gene$Ref <- str_trim(mutated_gene$Ref, side = "both")
mutated_gene$Alt <- str_trim(mutated_gene$Alt, side = "both")



#The followed part is based on one gene under different condition
#mutated gene information preprocess
#if mutation_position existed, get the mutated gene
#input the mutated information of gene from different conditons or strains

mutated_gene$complement_sign <- getSingleMatchParameter(gene_feature0$complement_sign,gene_feature0$locus_tag,mutated_gene$Gene2)
mutated_gene1 <- mutated_gene

for (i in seq(length(mutated_gene1$Chr))){
  if(mutated_gene1$complement_sign[i]){
  mutated_gene1$Ref[i] <- changeATCG(mutated_gene1$Ref[i])
  mutated_gene1$Alt[i] <- changeATCG(mutated_gene1$Alt[i])

  } else{
    mutated_gene1$Ref[i] <- mutated_gene1$Ref[i]
    mutated_gene1$Alt[i] <- mutated_gene1$Alt[i]
  }
}



#first run the program for each gene from different conditions or strains
#pre-process the gene annotation data before mutation mapping
#update the mutation information in the protein level
gene_list  <- unique(mutated_gene1$Gene2)
tt <- vector()
for (i in seq_along(gene_list)) {
  ss <- gene_list[i]
  gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature0)
  gene_snp[["pro_mutation_count"]] <- countMutationProtein(gene_name = ss, mutation_annotation = mutated_gene1, gene_snp0 = gene_snp)
  tt[i] <- sum(gene_snp[["pro_mutation_count"]])
  print(gene_snp)
}


#results analysis
num_gene_with_nsSNP <- tt[tt > 0]
num_nsSNP <- sum(num_gene_with_nsSNP)
protein_mutation <- data.frame(orf=gene_list,nsSNP=tt)



#-----------------------------------------------------------------------------------------
#analyze the mutation information with the structure
#-----------------------------------------------------------------------------------------
#first example
ss = 'YPR184W'
seq_from_3D <- 2:1534#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
dirForDistanceMatrix <- 'residue_distance/pdb_homo/2_1534_5d06.1.A_5b2453487f4bf94bf75ead43.pdb.txt'
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature0)
gene_snp[['pro_mutation_count']] <- countMutationProtein(gene_name = ss, mutation_annotation=mutated_gene1, gene_snp0 = gene_snp)
pos_mutation <- which(gene_snp[['pro_mutation_count']] != 0)

#residue_distance/pdb_homo/2_1534_5d06.1.A_5b2453487f4bf94bf75ead43.pdb.txt

#input the distance of all the pired residues
#ResidueDistance <- read_excel(dirForDistanceMatrix,col_names = FALSE) #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance <- read.table(dirForDistanceMatrix,sep = ",") #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance <- as.matrix(ResidueDistance)
dim(ResidueDistance)
ResidueDistance[1:5,1:5]
#obtain the mutation information for the structure
residueIn3D <- gene_snp[['protein']][seq_from_3D]
pos_mutation_3D <- gene_snp[['pro_mutation_count']][seq_from_3D]


#mutation position on structure and #mutation number on structure
pos_mutation_c <- which(pos_mutation_3D != 0)
seq0 <- 1:length(pos_mutation_3D) #seq0 is the coordinate of PDB structure
pos_count_num <- pos_mutation_3D[pos_mutation_c]


#calculate p_values using CLUPMS method
sample_standard1 <- sampleStand(pos_mutation_3D)
wap_original <- getTotalWAP(pos_mutation_c,sample_standard1,ResidueDistance)
wap_sample0 <- getSampleWAP(pos_mutation_c,sample_standard1,ResidueDistance, seq=seq0,n=10000)
plotNullDistribution(wap_sample0, wap_original0 = wap_original)
Strain_3D <- getPvalue(wap_original,wap_sample0)


# print the mutation information for the input SNP list contained in the protein 3D structure
pdbID <- '2_1534_5d06.1.A_5b2453487f4bf94bf75ead43'
SNP_list <- printSNPforGene(gene0 = 'YPR184W',
                SNPlist0 = mutated_gene1,
                gene_annotation0 = gene_feature0,
                pdbID0 = pdbID,
                sstart0 = 2,
                send0 = 1534)


#-------------------------------------------------------------------------------
#second example
#-------------------------------------------------------------------------------
ss = 'YMR246W'
seq_from_3D <- 39:691#this is the coordinated of original protein sequence and should changed into 3D structure coordinates
dirForDistanceMatrix <- 'residue_distance/pdb_homo/39_691_5mst.1.A_5b41c4d68fd6f9da68b53e00.pdb.txt'
gene_snp <- getGeneCoordinate(gene_name = ss, genesum = gene_feature0)
gene_snp[['pro_mutation_count']] <- countMutationProtein(gene_name = ss, mutation_annotation=mutated_gene1, gene_snp0 = gene_snp)
pos_mutation <- which(gene_snp[['pro_mutation_count']] != 0)


#input the distance of all the pired residues
#ResidueDistance <- read_excel(dirForDistanceMatrix,col_names = FALSE) #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance <- read.table(dirForDistanceMatrix,sep = ",") #in the followed calculation, the matrix dosen't have the col and row names
ResidueDistance <- as.matrix(ResidueDistance)


#obtain the mutation information for the structure
residueIn3D <- gene_snp[['protein']][seq_from_3D]
pos_mutation_3D <- gene_snp[['pro_mutation_count']][seq_from_3D]


#mutation position on structure and #mutation number on structure
pos_mutation_c <- which(pos_mutation_3D != 0)
seq0 <- 1:length(pos_mutation_3D) #seq0 is the coordinate of PDB structure
pos_count_num <- pos_mutation_3D[pos_mutation_c]


#calculate p_values using UPMS method
sample_standard1 <- sampleStand(pos_mutation_3D)
wap_original <- getTotalWAP(pos_mutation_c,sample_standard1,ResidueDistance)
wap_sample0 <- getSampleWAP(pos_mutation_c,sample_standard1,ResidueDistance, seq=seq0,n=10000)
plotNullDistribution(wap_sample0, wap_original0 = wap_original)
Strain_3D <- getPvalue(wap_original,wap_sample0)

