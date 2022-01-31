#' base = "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait/inst/extdata/genomic/Desulfotomaculum_kuznetsovii_dsm_6115_ASM21470v1_dna_toplevel"
#' x <- readFASTA(file.path(base, "translated.faa"))
#' x = x[1:3]
#' protcheck(x) # TRUE
#' protcheck(paste(x, "Z", sep = "")) # FALSE
protcheck <- function(x) {
  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
  )

  all(strsplit(x, split = "")[[1]] %in% AADict)
}

protdict_setup <- function() {
  protdict = list(AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"),
                  charged = c("D", "E", "K", "R"),
                  polar_uncharged = c("S","T","N","Q"),
                  thermolabile = c("H","Q","T"),
                  EK = c("E","K"),
                  QT = c("Q","T"),
                  QH = c("Q","H"),
                  EFMR = c("E","F","M","R"),
                  charged = c("D","E","K","R","H"),
                  polar = c("S","T","N","Q"),
                  hydrophobic = c("A","V","I","L","M","F","Y","W"),
                  percent_ERK = c("E","R","K"),
                  LK = c("L","K"),
                  Q = c("Q"),
                  IVYWREL = c("I","V","Y","W","R","E","L"),
                  IVWL = c("I","V","W","L"),
                  KVYWREP = c("K","V","W","R","E","P"),
                  GARP = c("G","A","R","P"),
                  MFILVWYERP = c("M","F","I","L","V","W","Y","E","R","P"),
                  ILVYER = c("I","L","V","Y","E","R"),
                  MILVWYER = c("M","I","L","V","W","Y","E","R"),
                  MFILVYERP = c("M","F","I","L","V","Y","E","R","P"),
                  FILVYERP = c("F","I","L","V","Y","E","R","P"),
                  ILVWYEHR = c("I","L","V","W","Y","E","H","R"),
                  FILVWYERP = c("F","I","L","V","W","Y","E","R","P"),
                  ILVWYGERKP = c("I","L","V","W","Y","G","E","R","K","P"),
                  ILVYEHR = c("I","L","V","Y","E","H","R")
                )
  return(protdict)
}



AA_counts <- function(prot_seqs) {
  AA_counts_matrix = do.call("rbind",lapply(prot_seqs, AA_count))
  AA_counts = apply(AA_counts_matrix, 2, sum)
  mean_length = sum(AA_counts)/length(prot_seqs)
  seq_count = length(prot_seqs)
  result <- list("AA_counts" = AA_counts, "mean_length" = mean_length, "seq_count" = seq_count)
  return(result)
}

AA_count <- function(prot_seq) {
  prot_seq = gsub("\\*", "", prot_seq)
  if (protcheck(prot_seq) == FALSE) {
    print(prot_seq, "\n")
    #stop("prot_seq has unrecognized amino acid type")
  }

  # 20 Amino Acid Abbrevation Dictionary from
  # https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties

  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
  )
  AAC <- summary(
    factor(strsplit(prot_seq, split = "")[[1]], levels = AADict)
  )
  return(AAC)
}

AA_freq <- function(AA_counts) {
  AA_freq = AA_counts/sum(AA_counts)
  return(AA_freq)
}

percent_charged <- function(AA_counts) {
  #calculate the fraction of the proteome that is charged
  AA_set = c("D", "E", "K", "R")
  total = sum(AA_counts)
  charged = sum(AA_counts[AA_set])
  return(charged/total)
}

percent_polar_uncharged <- function(AA_counts) {
  #calculate the fraction of the proteome that is polar-uncharged
  AA_set = c('S','T','N','Q')
  total = sum(AA_counts)
  polar = sum(AA_counts[AA_set])
  return(polar/total)
}

percent_thermolabile <- function(AA_counts) {
  #calculate the fraction of the proteome that is thermolabile AAs
  AA_set = c('H','Q','T')
  total = sum(AA_counts)
  HQT = sum(AA_counts[AA_set])
  return(HQT/total)
}

EKQH <- function(AA_counts) {
  #calculate the ratio of EK/QH in the proteome
  AA_set1 = c('E','K')
  AA_set2 = c('Q','H')
  ratio = sum(AA_counts[AA_set1])/sum(AA_counts[AA_set2])
  return(ratio)
}

EKQT <- function(AA_counts) {
  #calculate the ratio of EK/QT in the proteome
  AA_set1 = c('E','K')
  AA_set2 = c('Q','T')
  ratio = sum(AA_counts[AA_set1])/sum(AA_counts[AA_set2])
  return(ratio)
}

EFMR <- function(AA_counts) {
  #calculate the fraction of EFMR in the proteome
  AA_set = c('E','F','M','R')
  total = sum(AA_counts)
  EFMR = sum(AA_counts[AA_set])
  return(EFMR/total)
}

polar_charged <- function(AA_counts) {
  #calculate the ratio of polar-uncharged to polar-charged in the proteome
  AA_set1 = c('S','T','N','Q')
  AA_set2 = c('D','E','K','R','H')
  polar = sum(AA_counts[AA_set1])
  charged = sum(AA_counts[AA_set2])
  ratio = polar/charged
  return(ratio)
}

polar_hydrophobic <- function(AA_counts) {
  #calculate the ratio of polar-uncharged to polar-charged in the proteome
  AA_set1 = c('S','T','N','Q')
  AA_set2 = c('A','V','I','L','M','F','Y','W')
  polar = sum(AA_counts[AA_set1])
  hydrophobic = sum(AA_counts[AA_set2])
  ratio = polar/hydrophobic
  return(ratio)
}

percent_ERK <- function(AA_counts) {
  #calculate the fraction of the proteome that is thermolabile AAs
  AA_set = c('E','R','K')
  total = sum(AA_counts)
  ERK = sum(AA_counts[AA_set])
  return(ERK/total)
}

LKQ <- function(AA_counts) {
  #calculate the ratio of polar-uncharged to polar-charged in the proteome
  AA_set1 = c('L','K')
  AA_set2 = c('Q')
  LK = sum(AA_counts[AA_set1])
  Q = sum(AA_counts[AA_set2])
  ratio = LK/Q
  return(ratio)
}

percent_IVYWREL <- function(AA_counts) {
  #calculate the fraction of IVYWREL in the proteome
  AA_set = c('I','V','Y','W','R','E','L')
  total = sum(AA_counts)
  IVYWREL = sum(AA_counts[AA_set])
  return(IVYWREL/total)
}

percent_IVWL <- function(AA_counts) {
  #calculate the fraction of IVWL in the proteome
  AA_set = c('I','V','W','L')
  total = sum(AA_counts)
  IVWL = sum(AA_counts[AA_set])
  return(IVWL/total)
}

percent_KVYWREP <- function(AA_counts) {
  #calculate the fraction of KVYWREP in the proteome
  AA_set = c('K','V','W','R','E','P')
  total = sum(AA_counts)
  KVYWREP = sum(AA_counts[AA_set])
  return(KVYWREP/total)
}

percent_GARP <- function(AA_counts) {
  #calculate the fraction of GARP in the proteome
  AA_set = c('G','A','R','P')
  total = sum(AA_counts)
  GARP = sum(AA_counts[AA_set])
  return(GARP/total)
}

percent_MFILVWYERP <- function(AA_counts) {
  #calculate the fraction of MFILVWYERP in the proteome
  AA_set = c('M','F','I','L','V','W','Y','E','R','P')
  total = sum(AA_counts)
  MFILVWYERP = sum(AA_counts[AA_set])
  return(MFILVWYERP/total)
}

percent_ILVYER <- function(AA_counts) {
  #calculate the fraction of ILVYER in the proteome
  AA_set = c('I','L','V','Y','E','R')
  total = sum(AA_counts)
  ILVYER = sum(AA_counts[AA_set])
  return(ILVYER/total)
}

percent_MILVWYER <- function(AA_counts) {
  #calculate the fraction of MILVWYER in the proteome
  AA_set = c('M','I','L','V','W','Y','E','R')
  total = sum(AA_counts)
  MILVWYER = sum(AA_counts[AA_set])
  return(MILVWYER/total)
}

percent_MFILVYERP <- function(AA_counts) {
  #calculate the fraction of MILVWYER in the proteome
  AA_set = c('M','F','I','L','V','Y','E','R','P')
  total = sum(AA_counts)
  MFILVYERP = sum(AA_counts[AA_set])
  return(MFILVYERP/total)
}

percent_FILVYERP <- function(AA_counts) {
  #calculate the fraction of FILVYERP in the proteome
  AA_set = c('F','I','L','V','Y','E','R','P')
  total = sum(AA_counts)
  FILVYERP = sum(AA_counts[AA_set])
  return(FILVYERP/total)
}

percent_ILVWYEHR <- function(AA_counts) {
  #calculate the fraction of ILVWYEHR in the proteome
  AA_set = c('I','L','V','W','Y','E','H','R')
  total = sum(AA_counts)
  ILVWYEHR = sum(AA_counts[AA_set])
  return(ILVWYEHR/total)
}

percent_FILVWYERP <- function(AA_counts) {
  #calculate the fraction of FILVWYERP in the proteome
  AA_set = c('F','I','L','V','W','Y','E','R','P')
  total = sum(AA_counts)
  FILVWYERP = sum(AA_counts[AA_set])
  return(FILVWYERP/total)
}

percent_ILVWYGERKP <- function(AA_counts) {
  #calculate the fraction of ILVWYGERKP in the proteome
  AA_set = c('I','L','V','W','Y','G','E','R','K','P')
  total = sum(AA_counts)
  ILVWYGERKP = sum(AA_counts[AA_set])
  return(ILVWYGERKP/total)
}

percent_ILVYEHR <- function(AA_counts) {
  #calculate the fraction of ILVYEHR in the proteome
  AA_set = c('I','L','V','Y','E','H','R')
  total = sum(AA_counts)
  ILVYEHR = sum(AA_counts[AA_set])
  return(ILVYEHR/total)
}
