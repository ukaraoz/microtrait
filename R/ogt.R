
#file = "/Users/ukaraoz/Work/microtrait/code/github/microtrait/inst/extdata/2501651207_average_features.txt"
#features = read.table(file, check.names = F, sep = "\t", header = T)
#coefficients = models[[22]][5] %>% as.data.frame()
#setdiff(colnames(features), coefficients[,1])

#genome_file = "/Users/ukaraoz/Work/microtrait/code/github/microtrait/inst/extdata/genomic/2503283023.fna"
#result = run.prodigal(genome_file = genome_file, fa_file = tempfile(), faa_file = tempfile(), mode = "single")
##fastafile = "/opt/OGT_prediction-1.0.2/prediction/genomes/cyanobacterium_stanieri/2503283023.fa"
#cds_file = result$fa_file
#protein_file = result$faa_file
#kingdom = "bac"
#outdir = "/Users/ukaraoz/Work/microtrait/code/github/microtrait/inst/extdata/ogt.test"
extract_features <- function(genome_file, cds_file, protein_file, kingdom) {
  id = fs::path_file(genome_file)
  id = gsub("\\.fna$|\\.faa$|\\.fa$|", "", id, perl = T)

  # 1
  genome_genomicfeatures = analyze_genome(genome_file)
  write.table(rbind(c(Genome = id, genome_genomicfeatures)),
              file = file.path(outdir, paste0(id, "_genomic_features.txt")),
              row.names = F, col.names = T, sep = "\t", quote = F)
  # 2
  trna_results = run.tRNAscan(genome_file = genome_file)
  genome_tRNAfeatures = analyze_tRNA(trna_results$tRNA_fa)
  write.table(rbind(c(Genome = id, genome_tRNAfeatures)),
              file = file.path(outdir, paste0(id, "_tRNA_features.txt")),
              row.names = F, col.names = T, sep = "\t", quote = F)

  # 3
  rrna_results = run.barrnap(genome_file = genome_file, kingdom = kingdom)
  genome_rRNAfeatures = analyze_rRNA(rrna_results$rRNA_fa)
  write.table(rbind(c(Genome = id, genome_rRNAfeatures)),
              file = file.path(outdir, paste0(id, "_rRNA_features.txt")),
              row.names = F, col.names = T, sep = "\t", quote = F)
  # 4
  genomesize = genome_genomicfeatures["genomic Total Size"]
  genome_orffeatures = analyze_orfs(cds_file, genomesize)
  write.table(rbind(c(Genome = id, genome_orffeatures)),
              file = file.path(outdir, paste0(id, "_ORF_features.txt")),
              row.names = F, col.names = T, sep = "\t", quote = F)
  # 5
  genome_proteinfeatures = prot_analysis(protein_file)
  write.table(rbind(c(Genome = id, genome_proteinfeatures)),
              file = file.path(outdir, paste0(id, "_protein_features.txt")),
              row.names = F, col.names = T, sep = "\t", quote = F)

  features = data.frame(feature =
                         c(names(genome_genomicfeatures),
                          names(genome_tRNAfeatures),
                          names(genome_rRNAfeatures),
                          names(genome_orffeatures),
                          names(genome_proteinfeatures)),
                        value =
                          as.numeric(
                          c(genome_genomicfeatures,
                            genome_tRNAfeatures,
                            genome_rRNAfeatures,
                            genome_orffeatures,
                            genome_proteinfeatures))) %>%
      dplyr::as_tibble()

  model.select = models %>% dplyr::filter(rank == "superkingdom" & clade == "Bacteria")

  features %>%
    dplyr::left_join(model.select, by = c("feature" = "feature")) %>%
    dplyr::mutate(total = `value.x` * `value.y`) %>%
    dplyr::select(total) %>% sum()
    rowwise(feature) %>%
    summarise(m = sum(total))
  # setdiff(c(names(genome_genomicfeatures),
  #           names(genome_tRNAfeatures),
  #           names(genome_rRNAfeatures),
  #           names(genome_orffeatures),
  #           names(genome_proteinfeatures)),
  #         allcolnames)
  # write
}
