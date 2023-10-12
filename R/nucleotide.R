# genome_file = "/opt/OGT_prediction-1.0.2/prediction/genomes/cyanobacterium_stanieri/2503283023.fa"
# prodigal_file = "/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/mrna.fna"
# nucl_seq = readFASTA(prodigal_file)

#a = readDNAStringSet(prodigal_file)

#' @importFrom Biostrings readDNAStringSet oligonucleotideFrequency
analyze_genome <- function(fastafile) {
  results = list()
  dinucleotides = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

  genome = Biostrings::readDNAStringSet(fastafile)

  N_counts = Biostrings::oligonucleotideFrequency(genome, width = 1, as.prob = FALSE, simplify.as = "collapsed")
  N_freq = N_counts/sum(N_counts)
  names(N_freq) = paste0("genomic Nucleotide Fraction: ", names(N_freq))

  ###
  ## dinucleotide frequencies of the genome, dinucleotides with N not taken into account
  ###
  # this doesn't work here due to the way dinucleotides are counted
  # Biostrings::dinucleotideFrequency(genome, as.prob = TRUE)
  # suppress warning for odd length sequences
  dinucleotide_counts = suppressWarnings(lapply(as.character(genome), stringr::str_count, c(pattern = dinucleotides)))
  # apply(as.matrix(dinucleotide_counts), sum, 2); doesn't work here due to the returned class type
  dinucleotide_counts = do.call("rbind", dinucleotide_counts)
  colnames(dinucleotide_counts) = dinucleotides
  dinucleotide_counts = apply(dinucleotide_counts, 2, sum)
  dinucleotide_freqs = dinucleotide_counts/sum(dinucleotide_counts)
  names(dinucleotide_freqs) = paste0("genomic Dinucleotide Fraction: ", names(dinucleotide_freqs))

  ###
  ## percent GC of the genome
  ###
  GC = (N_counts["G"] + N_counts["C"])/sum(N_counts)
  names(GC) = "genomic GC"

  ###
  ## genome size
  ###
  total_size = sum(N_counts)
  names(total_size) = "genomic Total Size"

  ###
  ## J2 metric of the genome
  ###
  J2 = j2(genome)
  names(J2) = "genomic J2"

  result = c(N_freq,
             dinucleotide_freqs,
             GC,
             total_size,
             J2)
  result
}

#' @importFrom Biostrings readBStringSet letterFrequency
analyze_tRNA <- function(fastafile) {
  tRNA = Biostrings::readBStringSet(fastafile)

  results = list()
  N_counts = Biostrings::letterFrequency(tRNA, letters= c("A", "C", "G", "T"), collapse=T)

  N_freq = N_counts/sum(N_counts)
  GC = as.numeric((N_counts["G"] + N_counts["C"])/sum(N_counts))

  names(N_freq) = paste0("tRNA Nucleotide Fraction: ", names(N_freq))
  names(GC) = paste0("tRNA GC")
  result = c(N_freq, GC)
  result
}

#' @import Biostrings
analyze_rRNA <- function(fastafile) {
  # import rRNA
  # from Bio import SeqIO
  # data = SeqIO.index("/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/barrnap/barrnap_results_assigned_renamed.fa", "fasta")
  # results ={}
  # N_counts = rRNA.counter(data)
  # results['Nucleotide Fraction']=rRNA.nucleotide_freq(N_counts)
  # results['GC']=rRNA.GC(N_counts)
  # results2 = {}
  # for key in results.keys():
  #   if isinstance(results[key],dict):
  #     for subkey in results[key].keys():
  #       results2[key+': '+subkey]=results[key][subkey]
  # else:
  #   results2[key] = results[key]
  #fastafile = "/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/barrnap/barrnap_results_assigned_renamed.fa"
  rRNA = readDNAStringSet(fastafile)
  N_counts = oligonucleotideFrequency(rRNA, width = 1, as.prob = FALSE, simplify.as = "collapsed")
  #names(N_counts) = paste0("rRNA Nucleotide: ", names(N_counts))
  N_freq = N_counts/sum(N_counts)
  names(N_freq) = paste0("rRNA Nucleotide Fraction: ", names(N_freq))

  GC = as.numeric((N_counts["rRNA Nucleotide: G"] + N_counts["rRNA Nucleotide: C"])/sum(N_counts))
  names(GC) = "rRNA GC"

  result = c(N_freq, GC)
  result
}

#' @import Biostrings
analyze_orfs <- function(fastafile, genomesize) {
  # import ORFs
  # from Bio import SeqIO
  # data = SeqIO.index("/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/mrna.fna", "fasta")
  # results ={}
  # (N_counts,N_genes) = ORFs.counter(data)
  # results['Density'] = ORFs.density(N_genes,t_size)
  # results['Length'] = ORFs.length(N_counts,N_genes)
  # results['Coding_Noncoding'] = ORFs.coding_noncoding(N_counts,t_size)
  # results['Coding'] = ORFs.coding(N_counts,t_size)
  # results['Nucleotide Fraction'] = ORFs.nucleotide_freq(N_counts)
  # results['Codon'] = ORFs.codon(data)
  # results['Start Codon'] = ORFs.start_codon(data)
  # results['Stop Codon'] = ORFs.stop_codon(data)
  # results['GC'] = ORFs.GC(N_counts)
  # results['AG'] = ORFs.AG(N_counts)
  # results['Dinucleotide Fraction'] = ORFs.dinucleotide_freq(data)
  # results2 = {}
  # for key in results.keys():
  #   if isinstance(results[key],dict):
  #     for subkey in results[key].keys():
  #       results2[key+': '+subkey]=results[key][subkey]
  #   else:
  #     results2[key] = results[key]
  # genomesize = 3163381
  #fastafile = "/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/mrna.fna"
  orfs = Biostrings::readDNAStringSet(fastafile)
  N_counts = Biostrings::oligonucleotideFrequency(orfs, width = 1, as.prob = FALSE, simplify.as = "collapsed")
  N_genes = length(orfs)
  density = N_genes/genomesize

  # calculate average ORF length
  length = sum(N_counts)/N_genes

  # calculate ratio of coding/noncoding
  coding_noncoding = sum(N_counts)/((2*genomesize) - sum(N_counts))

  # calculate fraction of the genome that is coding
  coding = sum(N_counts)/genomesize

  # calculate the nucleotide fraction
  nucleotide_freq = N_counts/sum(N_counts)

  # calculate the fraction of each codon in the ORFs
  codons = c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
             "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
             "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
             "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT")
  codon_counts = apply(Biostrings::oligonucleotideFrequency(orfs, width = 3, step = 3), 2, sum)
  codon_counts = codon_counts/sum(codon_counts)

  # calculate the fraction of each start codon
  orfs_startcodon = Biostrings::subseq(orfs, start = 1, end = 3)
  startcodon = apply(Biostrings::oligonucleotideFrequency(orfs_startcodon, width = 3), 2, sum)
  startcodon = startcodon[c("ATG","GTG","TTG")]
  startcodon_freq = startcodon/sum(startcodon)

  # calculate the fraction of each stop codon
  temp = as.character(orfs)
  orfs_stopcodon = Biostrings::DNAStringSet(substr(temp, nchar(temp) - 2, nchar(temp)))
  stopcodon = apply(Biostrings::oligonucleotideFrequency(orfs_stopcodon, width = 3), 2, sum)
  stopcodon = stopcodon[c("TAA","TAG","TGA")]
  stopcodon_freq = stopcodon/sum(stopcodon)

  GC = as.numeric((N_counts["G"] + N_counts["C"])/sum(N_counts))
  AG = as.numeric((N_counts["G"] + N_counts["A"])/sum(N_counts))

  dinucleotides = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
  dinucleotide_counts = suppressWarnings(lapply(as.character(orfs), stringr::str_count, c(pattern = dinucleotides)))
  dinucleotide_counts = do.call("rbind", dinucleotide_counts)
  colnames(dinucleotide_counts) = dinucleotides
  dinucleotide_counts = apply(dinucleotide_counts, 2, sum)
  dinucleotide_freq = dinucleotide_counts/sum(dinucleotide_counts)
  rownames(dinucleotide_freq) = NULL

  names(density) = "ORF Density"
  names(length) = "ORF Length"
  names(coding_noncoding) = "ORF Coding_Noncoding"
  names(coding) = "ORF Coding"
  names(nucleotide_freq) = paste0("ORF Nucleotide Fraction: ", names(nucleotide_freq))
  names(codon_counts) = paste0("ORF Codon: ", names(codon_counts))
  names(startcodon_freq) = paste0("ORF Start Codon: ", names(startcodon_freq))
  names(stopcodon_freq) = paste0("ORF Stop Codon: ", names(stopcodon_freq))
  names(GC) = "ORF GC"
  names(AG) = "ORF AG"
  names(dinucleotide_freq) = paste0("ORF Dinucleotide Fraction: ", names(dinucleotide_freq))

  result = c(density, length, coding_noncoding, coding,
             nucleotide_freq, codon_counts, startcodon_freq, stopcodon_freq,
             GC, AG, dinucleotide_freq)
  return(result)
}

j2 <- function(genome) {
  # x is of type XStringSet
  # calculate the J2 metric of the genome
  dinucleotides = c("TT","CC","TC","CT","AA","GG","AG","GA",
                    "TA","TG","CA","CG", "AT","AC","GT","GC")
  dinucleotide_counts = suppressWarnings(lapply(as.character(genome), stringr::str_count, c(pattern = dinucleotides)))
  dinucleotide_counts = do.call("rbind", dinucleotide_counts)
  colnames(dinucleotide_counts) = dinucleotides
  dinucleotide_counts = apply(dinucleotide_counts, 2, sum)

  YY = sum(dinucleotide_counts[c("TT","CC","TC","CT")])
  RR = sum(dinucleotide_counts[c("AA","GG","AG","GA")])
  YR = sum(dinucleotide_counts[c("TA","TG","CA","CG")])
  RY = sum(dinucleotide_counts[c("AT","AC","GT","GC")])
  total = YY+RR+YR+RY
  YY = YY/total
  RR = RR/total
  YR = YR/total
  RY = RY/total
  return(YY+RR-YR-RY)
}
