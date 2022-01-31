# genome_file = "/opt/OGT_prediction-1.0.2/prediction/genomes/cyanobacterium_stanieri/2503283023.fa"
# prodigal_file = "/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/mrna.fna"
# nucl_seq = readFASTA(prodigal_file)

#a = readDNAStringSet(prodigal_file)
analyze_genome <- function(fastafile) {
  # import genomic
  # from Bio import SeqIO
  # data = SeqIO.index("/opt/OGT_prediction-1.0.2/prediction/genomes/cyanobacterium_stanieri/2503283023.fa", "fasta")
  # #data = SeqIO.index("/opt/OGT_prediction-1.0.2/prediction/genomes/cyanobacterium_stanieri/test.fa", "fasta")
  # results = {}
  # N_counts = genomic.counter(data)
  # results['Nucleotide Fraction']=genomic.nucleotide_freq(N_counts)
  # results['Dinucleotide Fraction']=genomic.dinucleotide_freq(data)
  # results['GC']=genomic.GC(N_counts)
  # results['Total Size']=genomic.t_size(N_counts)
  # results['J2']=genomic.j2(data)
  # results2 = {}
  # for key in results.keys():
  #   if isinstance(results[key],dict):
  #     for subkey in results[key].keys():
  #       results2[key+': '+subkey]=results[key][subkey]
  #   else:
  #     results2[key] = results[key]
  #fastafile = "/opt/OGT_prediction-1.0.2/prediction/genomes/cyanobacterium_stanieri/2503283023.fa"
  dinucleotides = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
  genome = Biostrings::readDNAStringSet(fastafile)
  results = list()
  N_counts = Biostrings::oligonucleotideFrequency(genome, width = 1, as.prob = FALSE, simplify.as = "collapsed")
  N_freq = N_counts/sum(N_counts)
  names(N_freq) = paste0("genomic Nucleotide Fraction: ", names(N_freq))

  # this doesn't generate same results due to the way dinucleotides are counted
  # Biostrings::dinucleotideFrequency(genome, as.prob = TRUE)
  dinucleotide_counts = stringr::str_count(as.character(genome), dinucleotides)
  names(dinucleotide_counts) = dinucleotides
  dinucleotide_freqs = dinucleotide_counts/sum(dinucleotide_counts)
  names(dinucleotide_freqs) = paste0("genomic Dinucleotide Fraction: ", names(dinucleotide_freqs))

  GC = (N_counts["G"] + N_counts["C"])/sum(N_counts)
  names(GC) = "genomic GC"

  total_size = sum(N_counts)
  names(total_size) = "genomic Total Size"

  J2 = j2(genome)
  names(J2) = "genomic J2"

  result = c(N_freq,
             dinucleotide_freqs,
             GC,
             total_size,
             J2)
  result
}

analyze_tRNA <- function(fastafile) {
  # import tRNA
  # from Bio import SeqIO
  # data = SeqIO.index("/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/trnascan_result.fa", "fasta")
  # results ={}
  # N_counts = tRNA.counter(data)
  # results['Nucleotide Fraction']=tRNA.nucleotide_freq(N_counts)
  # results['GC']=tRNA.GC(N_counts)
  # results2 = {}
  # for key in results.keys():
  #   if isinstance(results[key],dict):
  #     for subkey in results[key].keys():
  #       results2[key+': '+subkey]=results[key][subkey]
  # else:
  #   results2[key] = results[key]
  #fastafile = "/opt/OGT_prediction-1.0.2/prediction/output/genomes/cyanobacterium_stanieri/2503283023/trnascan_result_masked.fa"
  tRNA = Biostrings::readDNAStringSet(fastafile)
  results = list()
  N_counts = Biostrings::oligonucleotideFrequency(tRNA, width = 1, as.prob = FALSE, simplify.as = "collapsed")
  N_freq = N_counts/sum(N_counts)
  names(N_freq) = paste0("tRNA Nucleotide Fraction: ", names(N_freq))

  GC = as.numeric((N_counts["tRNA Nucleotide Fraction: G"] + N_counts["tRNA Nucleotide Fraction: C"])/sum(N_counts))
  names(GC) = paste0("tRNA GC")

  result = c(N_freq, GC)
  result
}

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
  orfs = readDNAStringSet(fastafile)
  N_counts = oligonucleotideFrequency(orfs, width = 1, as.prob = FALSE, simplify.as = "collapsed")
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
  codon_counts = apply(oligonucleotideFrequency(orfs, width = 3, step = 3), 2, sum)
  codon_counts = codon_counts/sum(codon_counts)

  # calculate the fraction of each start codon
  orfs_startcodon = subseq(orfs, start = 1, end = 3)
  startcodon = apply(oligonucleotideFrequency(orfs_startcodon, width = 3), 2, sum)
  startcodon = startcodon[c("ATG","GTG","TTG")]
  startcodon_freq = startcodon/sum(startcodon)

  # calculate the fraction of each stop codon
  temp = as.character(orfs)
  orfs_stopcodon = DNAStringSet(substr(temp, nchar(temp) - 2, nchar(temp)))
  stopcodon = apply(oligonucleotideFrequency(orfs_stopcodon, width = 3), 2, sum)
  stopcodon = stopcodon[c("TAA","TAG","TGA")]
  stopcodon_freq = stopcodon/sum(stopcodon)

  GC = as.numeric((N_counts["G"] + N_counts["C"])/sum(N_counts))
  AG = as.numeric((N_counts["G"] + N_counts["A"])/sum(N_counts))

  dinucleotides = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
  dinucleotide_counts = lapply(orfs, function(x) {stringr::str_count(as.character(x), dinucleotides)})
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

j2 <- function(x) {
  # x is of type XStringSet
  # calculate the J2 metric of the genome
  YY = sum(stringr::str_count(as.character(x), c("TT","CC","TC","CT")))
  RR = sum(stringr::str_count(as.character(x), c("AA","GG","AG","GA")))
  YR = sum(stringr::str_count(as.character(x), c("TA","TG","CA","CG")))
  RY = sum(stringr::str_count(as.character(x), c("AT","AC","GT","GC")))
  total = YY+RR+YR+RY
  YY = YY/total
  RR = RR/total
  YR = YR/total
  RY = RY/total
  return(YY+RR-YR-RY)
}

