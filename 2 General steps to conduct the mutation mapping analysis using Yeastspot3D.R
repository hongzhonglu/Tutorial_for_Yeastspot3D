# Note:
# This script is showing how to conduct the mutation mapping analysis for the small SNP dataset
data('gene_feature0')
data('snp_data')
mutated_gene <- annotateSNP(snp_input = snp_data, gene_feature = gene_feature0)


#-------------------------------------------------
# Mutation enrichment analysis
#------------------------------------------------
# first example
data('ResidueDistance_YPR184W')
mutated_gene1 <- filter(mutated_gene, Gene2 == 'YPR184W')
result0 <- clumpsAnalysis(gene0 = 'YPR184W',
                          SNPlist0 = mutated_gene1,
                          gene_annotation0 = gene_feature0,
                          pdb = ResidueDistance_YPR184W,
                          sstart0 = 2,
                          send0 = 1534,
                          input_dir= FALSE)

# print the mutation information for the input SNP list contained in the protein 3D structure
pdbID <- '2_1534_5d06.1.A_5b2453487f4bf94bf75ead43'
SNP_list <- printSNPforGene(gene0 = 'YPR184W',
                SNPlist0 = mutated_gene1,
                gene_annotation0 = gene_feature0,
                pdbID0 = pdbID,
                sstart0 = 2,
                send0 = 1534)


# second example
data('ResidueDistance_YMR246W')
mutated_gene1 <- filter(mutated_gene, Gene2 == 'YMR246W')
result0 <- clumpsAnalysis(gene0 = 'YMR246W',
                          SNPlist0 = mutated_gene1,
                          gene_annotation0 = gene_feature0,
                          pdb = ResidueDistance_YMR246W,
                          sstart0 = 39,
                          send0 = 691,
                          input_dir= FALSE)

pdbID <- '39_691_5mst.1.A_5b41c4d68fd6f9da68b53e00'
SNP_list <- printSNPforGene(gene0 = 'YMR246W',
                            SNPlist0 = mutated_gene1,
                            gene_annotation0 = gene_feature0,
                            pdbID0 = pdbID,
                            sstart0 = 39,
                            send0 = 691)


#-------------------------------------------------
# Mutation hot spot analysis
#------------------------------------------------
# run the function
data('snp_YBR046C')
data('ResidueDistance_YBR046C')
outfile0 <- 'result/hot_spot_analysis'
dir.create(outfile0)

hotSpotAnalysis(
  gene0 = "YBR046C",
  SNPlist0 = snp_YBR046C,
  gene_annotation0 = gene_feature0,
  pdb = ResidueDistance_YBR046C,
  sstart0 = 5, # coordinate of orginal protein residues sequence
  send0 = 333,     # coordinate of orginal protein residues sequence
  qstart0 =1 , # coordinate of protein residues sequence in pdb file
  qend0 = 329,     # coordinate of protein residues sequence in pdb file
  result_dir = outfile0,
  input_dir=FALSE
)
