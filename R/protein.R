protdict_setup <- function() {
  protdict = list(AA = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"),
                  charged = c("D", "E", "K", "R"),
                  chargedwH = c("D","E","K","R","H"),
                  polaruncharged = c("S","T","N","Q"),
                  thermolabile = c("H","Q","T"),
                  EK = c("E","K"),
                  QT = c("Q","T"),
                  QH = c("Q","H"),
                  EFMR = c("E","F","M","R"),
                  ERK = c("E","R","K"),
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

prot_analysis <- function(fastafile) {
  #fastafile = "/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/translated.faa"
  #import protein
  #import pprint
  #from Bio import SeqIO
  #protein_seqs = SeqIO.index("/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/translated.faa", "fasta")
  # AA_counts = protein.AA_count(protein_seqs)
  #result = protein.analysis(protein_seqs)
  #pp = pprint.PrettyPrinter(indent=4)
  #pp.pprint(result)

  protdict = protdict_setup()
  proteins = readAAStringSet(fastafile)
  AA_count_matrix = letterFrequency(proteins, letters = protdict$AA)
  AA_counts = apply(AA_count_matrix, 2, sum)

  mean_length = sum(AA_counts)/length(proteins)
  names(mean_length) = paste0("protein Length", names(mean_length))

  AA_freq = AA_counts/(sum(AA_counts))
  names(AA_freq) = paste0("protein AA: ", names(AA_freq))

  ratio_charged = sum(AA_counts[protdict$charged]) / sum(AA_counts)
  names(ratio_charged) = paste0("protein ", "Charged")

  ratio_polaruncharged = sum(AA_counts[protdict$polaruncharged]) / sum(AA_counts)
  names(ratio_polaruncharged) = paste0("protein ", "Polar-Uncharged")

  ratio_thermolabile = sum(AA_counts[protdict$thermolabile]) / sum(AA_counts)
  names(ratio_thermolabile) = paste0("protein ", "Thermolabile")

  EK_over_QH = sum(AA_counts[protdict$EK]) / sum(AA_counts[protdict$QH])
  names(EK_over_QH) = paste0("protein ", "EK_QH")

  EK_over_QT = sum(AA_counts[protdict$EK]) / sum(AA_counts[protdict$QT])
  names(EK_over_QT) = paste0("protein ", "EK_QT")

  ratio_EFMR = sum(AA_counts[protdict$EFMR]) / sum(AA_counts)
  names(ratio_EFMR) = paste0("protein ", "EFMR")

  polar_over_charged = sum(AA_counts[protdict$polar]) / sum(AA_counts[protdict$chargedwH])
  names(polar_over_charged) = paste0("protein ", "Polar_Charged")

  polar_over_hydrophobic = sum(AA_counts[protdict$polar]) / sum(AA_counts[protdict$hydrophobic])
  names(polar_over_hydrophobic) = paste0("protein ", "Polar_Hydrophobic")

  ratio_ERK = sum(AA_counts[protdict$ERK]) / sum(AA_counts)
  names(ratio_ERK) = paste0("protein ", "ERK")

  LK_over_Q = sum(AA_counts[protdict$LK]) / sum(AA_counts[protdict$Q])
  names(LK_over_Q) = paste0("protein ", "LK_Q")

  ratio_IVYWREL = sum(AA_counts[protdict$IVYWREL]) / sum(AA_counts)
  names(ratio_IVYWREL) = paste0("protein ", "IVYWREL")

  ratio_IVWL = sum(AA_counts[protdict$IVWL]) / sum(AA_counts)
  names(ratio_IVWL) = paste0("protein ", "IVWL")

  ratio_KVYWREP = sum(AA_counts[protdict$KVYWREP]) / sum(AA_counts)
  names(ratio_KVYWREP) = paste0("protein ", "KVYWREP")

  ratio_GARP = sum(AA_counts[protdict$GARP]) / sum(AA_counts)
  names(ratio_GARP) = paste0("protein ", "GARP")

  ratio_MFILVWYERP = sum(AA_counts[protdict$MFILVWYERP]) / sum(AA_counts)
  names(ratio_MFILVWYERP) = paste0("protein ", "MFILVWYERP")

  ratio_ILVYER = sum(AA_counts[protdict$ILVYER]) / sum(AA_counts)
  names(ratio_ILVYER) = paste0("protein ", "ILVYER")

  ratio_MILVWYER = sum(AA_counts[protdict$MILVWYER]) / sum(AA_counts)
  names(ratio_MILVWYER) = paste0("protein ", "MILVWYER")

  ratio_MFILVYERP = sum(AA_counts[protdict$MFILVYERP]) / sum(AA_counts)
  names(ratio_MFILVYERP) = paste0("protein ", "MFILVYERP")

  ratio_FILVYERP = sum(AA_counts[protdict$FILVYERP]) / sum(AA_counts)
  names(ratio_FILVYERP) = paste0("protein ", "FILVYERP")

  ratio_ILVWYEHR = sum(AA_counts[protdict$ILVWYEHR]) / sum(AA_counts)
  names(ratio_ILVWYEHR) = paste0("protein ", "ILVWYEHR")

  ratio_FILVWYERP = sum(AA_counts[protdict$FILVWYERP]) / sum(AA_counts)
  names(ratio_FILVWYERP) = paste0("protein ", "FILVWYERP")

  ratio_ILVWYGERKP = sum(AA_counts[protdict$ILVWYGERKP]) / sum(AA_counts)
  names(ratio_ILVWYGERKP) = paste0("protein ", "ILVWYGERKP")

  ratio_ILVYEHR = sum(AA_counts[protdict$ILVYEHR]) / sum(AA_counts)
  names(ratio_ILVYEHR) = paste0("protein ", "ILVYEHR")

  result = c(mean_length,
             AA_freq,
             ratio_charged,
             ratio_polaruncharged,
             ratio_thermolabile,
             EK_over_QH,
             EK_over_QT,
             ratio_EFMR,
             polar_over_charged,
             polar_over_hydrophobic,
             ratio_ERK,
             LK_over_Q,
             ratio_IVYWREL,
             ratio_IVWL,
             ratio_KVYWREP,
             ratio_GARP,
             ratio_MFILVWYERP,
             ratio_ILVYER,
             ratio_MILVWYER,
             ratio_MFILVYERP,
             ratio_FILVYERP,
             ratio_ILVWYEHR,
             ratio_FILVWYERP,
             ratio_ILVWYGERKP,
             ratio_ILVYEHR)
  result
}

