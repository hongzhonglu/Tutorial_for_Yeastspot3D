#-------------------------------------------------------------
# this code is used to extract the resolution from pdb files
#-------------------------------------------------------------
library(stringr)
library(tidyverse)
protein.ResolutionEX <- function(pdbdir) {
  # input a pdb file
  # output the resolution for the pdb file
  # experiment pdb file
  pdb <- scan(pdbdir, sep = "\n", what = "complex")
  ss1 <- which(str_detect(pdb, "REMARK   2") == TRUE)
  pdb1 <- pdb[ss1]
  ss2 <- which(str_detect(pdb1, "RESOLUTION") == TRUE)
  pdb2 <- pdb1[ss2]
  pdb3 <- str_split(pdb2, "RESOLUTION")
  # obtain the resolution
  resolution <- pdb3[[1]][2] %>%
    str_extract(.,"\\d+\\.*\\d*")
  return(resolution)
}

#test
infile <- "data/"
pdbid <- "2nrn.pdb"
pdbdir <- paste(infile, pdbid, sep = "")
protein.ResolutionEX(pdbdir)


pdbid <- "2lhh.pdb"
pdbdir <- paste(infile, pdbid, sep = "")
protein.ResolutionEX(pdbdir)


#----------------------------------------------------------------------------
# this code is used to extract the residue sequence from experiment pdb files
#----------------------------------------------------------------------------
library(bio3d)
library(seqinr)
infile <- "data/"
pdbid <- '6cp6.pdb'

pdbdir <- paste(infile, pdbid, sep = "")

protein.SequenceEX <- function(pdbdir) {
  pdb <- read.pdb(pdbdir)
  atom1 <- pdb$atom
  chainAll <- unique(atom1$chain)
  seq_list <- list()
  for (j in chainAll) {
    cat("chainID")
    print(j)
    pdbid0 <- str_replace_all(pdbid, "pdb", j)
    # should remove the legend not belong to amino acids
    atom_choose <- atom1[atom1$chain == j & atom1$type == "ATOM", ]
    atom_choose_o <- select(atom_choose, resid, chain, resno)
    for (i in seq_along(atom_choose$resid)) {
      atom_choose$resid[i] <- str_to_title(atom_choose$resid[i])
      if (is.na(a(atom_choose$resid[i])) & atom_choose$resid[i] == "Fme") {
        # sometimes a molecular can not be a amino acid residue
        atom_choose_o$resid[i] <- "M"
      } else if (is.na(a(atom_choose$resid[i])) & atom_choose$resid[i] != "Fme") {
        atom_choose_o$resid[i] <- "*"
      } else {
        # change three letter into one
        atom_choose_o$resid[i] <- a(atom_choose$resid[i])
      }
    }
    atom_choose_o$not_duplicated <- !duplicated(atom_choose_o$resno)
    seq0 <- atom_choose_o[atom_choose_o$not_duplicated == TRUE, ] %>%
      select(., resid)
    seq0 <- paste(seq0$resid, collapse = "")
    print(seq0)
    seq_list[[pdbid0]] <- seq0
  }
  
  return(seq_list)
}

protein.SequenceEX(pdbdir)


#----------------------------------------------------------------------------
# this code is used to blast analysis of amino acids
#----------------------------------------------------------------------------
library(metablastr)
# run blastn (nucleotide to nucleotide search) between example query and subject sequences
# fistly download the blast and installed it
# wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.2.28/ncbi-blast-2.2.28+.dmg
blast_test <- blast_nucleotide_to_nucleotide(
  query   = system.file('seqs/qry_nn.fa', package = 'metablastr'),
  subject = system.file('seqs/sbj_nn.fa', package = 'metablastr'),
  output.path = tempdir(),
  db.import  = FALSE)
# look at BLAST results
blast_test



#----------------------------------------------------------------------------
# this code is used to calculate the distance
#----------------------------------------------------------------------------
infile <- "data/"
pdbid <- '6cp6.pdb'
pdbdir <- paste(infile, pdbid, sep = "")

 
protein.ResidueDistance <- function(pdbdir,chainID) {
  pdb <- read.pdb(pdbdir)
  sele <- atom.select(pdb, type = "ATOM", "calpha", chain = chainID)
  pdb0 <- trim.pdb(pdb, sele)
  newPDB <- pdb0$atom
  k <- dm(pdb0, inds = "calpha")
  filled.contour(k, nlevels = 10)
  return(k)
}

protein.ResidueDistance(pdbdir, chainID = 'K')




  