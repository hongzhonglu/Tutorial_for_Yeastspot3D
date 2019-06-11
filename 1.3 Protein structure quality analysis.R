# This code is mainly used for the protein structure quality analysis
# will refine from 2019-1-28



#-----------------------------------------------------------------
# step 1 read data
#-----------------------------------------------------------------
# This file should have at least four parameters for each homology protein structures
# These four parameters are
# qmean
# Seq_similarity
# Seq-identity
# Resolution
pdb_homo_test1 <- read.table('data/pdb_HOMO_test.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)

pdbHomoFilter <- function(pdb_homo, qmean0=-4, identity0=25, similarity0=0.31, resolution0=3.4 ){
  colname0 <- colnames(pdb_homo)
  #check whether the four parameters as the column names
  key_para_filter <- c('qmean','Seq_Identity','Seq_similarity','Resolution')
  if (all(key_para_filter %in% colname0)){
    print('all four key parameters exist')
    print('Conduct the filteration for homology PDB files')
    pdb_homo$qmean <- as.numeric(pdb_homo$qmean)
    pdb_homo$Seq_Identity <- as.numeric(pdb_homo$Seq_Identity)
    pdb_homo$Seq_similarity <- as.numeric(pdb_homo$Seq_similarity)
    pdb_homo$Resolution <- as.numeric(pdb_homo$Resolution)
    pdb_homo_filter <- filter(pdb_homo, pdb_homo$qmean >= qmean0 &
                                pdb_homo$Seq_Identity >= identity0 &
                                pdb_homo$Seq_similarity >=similarity0 &
                                pdb_homo$Resolution <= resolution0 )
  } else{
    print('Please prepare the detailed parameters for the quality analysis')
  }
}

pdb_homo_filter1 <- pdbHomoFilter(pdb_homo=pdb_homo_test1)



# This file should have at least four parameters for each homology protein structures
# These four parameters are
# pident
# mismatch
# Resolution
pdb_EX_test <- read.table('data/pdb_EX_test.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pdb_EX_test$pident <-  sample(95:100, 3365, replace=T)
pdb_EX_test$mismatch <-  sample(0:1, 3365, replace=T)

pdbExFilter <- function(pdb_EX, pident0=100, mismatch0=0, resolution0=3.4 ){
  colname0 <- colnames(pdb_EX)
  #check whether the four parameters as the column names
  key_para_filter <- c('pident','mismatch','Resolution')
  if (all(key_para_filter %in% colname0)){
    print('all three key parameters exist')
    print('Conduct the filteration for experimental PDB files')
    pdb_EX$pident <- as.numeric(pdb_EX$pident)
    pdb_EX$mismatch <- as.numeric(pdb_EX$mismatch)
    pdb_EX$Resolution <- as.numeric(pdb_EX$Resolution)
    pdb_EX_filter <- filter(pdb_EX, pdb_EX$pident >= pident0 &
                              pdb_EX$mismatch <= mismatch0 &
                              pdb_EX$Resolution <= resolution0 )
  } else{
    print('Please prepare the detailed parameters for the quality analysis')
  }
}
pdb_EX_filter1 <- pdbExFilter(pdb_EX = pdb_EX_test)








