# Note:
# This script is showing how to conduct the mutation mapping analysis for the small SNP dataset

# SNP data preparation
dir0 <- "data/snp_adaption_to_high_ethanol.XLS"
snp_data <- read_excel(dir0)
mutated_gene1 <- annotateSNP(snp_input = snp_data, gene_feature = gene_feature0)


#-------------------------------------------------
# Mutation enrichment analysis
#------------------------------------------------
# first example
ss0 <- 'YPR184W'
p1 <- 2
p2 <- 1534
distance_dir <- 'residue_distance/pdb_homo/2_1534_5d06.1.A_5b2453487f4bf94bf75ead43.pdb.txt'
result0 <- clumpsAnalysis(gene0 = ss0,
                          SNPlist0 = mutated_gene1,
                          gene_annotation0 = gene_feature0,
                          pdb_dir = distance_dir,
                          sstart0 = p1,
                          send0 = p2)



# print the mutation information for the input SNP list contained in the protein 3D structure
pdbID <- '2_1534_5d06.1.A_5b2453487f4bf94bf75ead43'
SNP_list <- printSNPforGene(gene0 = 'YPR184W',
                SNPlist0 = mutated_gene1,
                gene_annotation0 = gene_feature0,
                pdbID0 = pdbID,
                sstart0 = 2,
                send0 = 1534)

# second example
ss0 <- 'YMR246W'
p1 <- 39
p2 <- 691
distance_dir <- 'residue_distance/pdb_homo/39_691_5mst.1.A_5b41c4d68fd6f9da68b53e00.pdb.txt'
result0 <- clumpsAnalysis(gene0 = ss0,
                          SNPlist0 = mutated_gene1,
                          gene_annotation0 = gene_feature0,
                          pdb_dir = distance_dir,
                          sstart0 = p1,
                          send0 = p2)


#-------------------------------------------------
# Mutation hot spot analysis
#------------------------------------------------
ss0 <- 'YPR184W'
p1 <- 2
p2 <- 1534
distance_dir <- 'residue_distance/pdb_homo/2_1534_5d06.1.A_5b2453487f4bf94bf75ead43.pdb.txt'
outfile0 <- 'result/hot_spot_analysis'
dir.create(outfile0)
hotSpotAnalysis(
  gene0 = ss0,
  SNPlist0 = mutated_gene1,
  gene_annotation0 = gene_feature0,
  pdb_dir = distance_dir,
  sstart0 = p1,
  send0 = p2,
  result_dir = outfile0)
