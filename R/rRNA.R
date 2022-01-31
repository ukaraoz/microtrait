extract_rRNA <- function(genome_file, kingdom, threads = 1) {
  #Name=16S_rRNA;product=16S ribosomal RNA
  #genome_file = "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait/inst/extdata/genomic/2503283023.fna"
  barrnap_result = run.barrnap(genome_file, kingdom = kingdom, threads = threads)
  bed16S_file = barrnap_result$tempfile_16S
  fa16S_file = run.bedtools(genome_file, bed16S_file = bed16S_file)
  return(fa16S_file)
}

#fa16S_file = extract_rRNA(genome_file, kingdom = "bac")
#nucl_analysis(fa16S_file)

#' base = "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait/inst/extdata/genomic/Desulfotomaculum_kuznetsovii_dsm_6115_ASM21470v1_dna_toplevel"
#' x <- readFASTA(file.path(base, "translated.faa"))
#' x = x[1:3]
#' protcheck(x) # TRUE
#' protcheck(paste(x, "Z", sep = "")) # FALSE
nuclcheck <- function(x) {
  NuclDict <- c("A", "G", "T", "C")

  all(strsplit(x, split = "")[[1]] %in% NuclDict)
}

nucl_analysis <- function(fastafile) {
  #fastafile = "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait/inst/extdata/genomic/Desulfotomaculum_kuznetsovii_dsm_6115_ASM21470v1_dna_toplevel/translated.faa"
  seqs <- readFASTA(fastafile)
  temp = nucl_counts(seqs)
  nucl_counts = temp$nucl_counts
  mean_length = temp$mean_length
  seq_count = temp$seq_count

  nucl_freq = nucl_freq(nucl_counts)
  GC = GC_count(nucl_counts)

  result = list(seq_count = seq_count,
                mean_length = mean_length,
                nucl_freq = nucl_freq,
                GC = GC
                )
  result
}

nucl_counts <- function(nucl_seqs) {
  nucl_counts_matrix = do.call("rbind",lapply(nucl_seqs, nucl_count))
  nucl_counts = apply(nucl_counts_matrix, 2, sum)
  mean_length = sum(nucl_counts)/length(nucl_seqs)
  seq_count = length(nucl_seqs)
  result <- list("nucl_counts" = nucl_counts, "mean_length" = mean_length, "seq_count" = seq_count)
  return(result)
}

nucl_count <- function(nucl_seq) {
  nucl_seq = gsub("\\*", "", nucl_seq)
  if (nuclcheck(nucl_seq) == FALSE) {
    print(nucl_seq, "\n")
    #stop("nucl_seq has unrecognized amino acid type")
  }

  NuclDict <- c("A", "G", "T", "C")
  NuclC <- summary(
    factor(strsplit(nucl_seq, split = "")[[1]], levels = NuclDict)
  )
  return(NuclC)
}

GC_count <- function(nucl_counts) {
  #calculate the GC fraction
  nucl_set = c('G','C')
  total = sum(nucl_counts)
  GC = sum(nucl_counts[nucl_set])
  return(GC/total)
}


