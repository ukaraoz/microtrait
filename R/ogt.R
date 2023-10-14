#' @param genome_file genome_file
#' @param cds_file cds_file
#' @param proteins_file proteins_file
#' @importFrom Biostrings readDNAStringSet
#' @importFrom dplyr pull as_tibble left_join mutate select filter
#' @import gRodon
#' @returns features
#' @export extract_features
extract_features <- function(genome_file, cds_file, proteins_file, tRNA = TRUE) {
  id = fs::path_file(genome_file)
  id = gsub("\\.fna$|\\.faa$|\\.fa$|", "", id, perl = T)

  # 1
  message("Extracting genome level features")
  genome_genomicfeatures = analyze_genome(genome_file)
  # write.table(rbind(c(Genome = id, genome_genomicfeatures)),
  #             file = file.path(outdir, paste0(id, "_genomic_features.txt")),
  #             row.names = F, col.names = T, sep = "\t", quote = F)
  # 3
  #rrna_results = run.barrnap(genome_file = genome_file, kingdom = kingdom)
  #genome_rRNAfeatures = analyze_rRNA(rrna_results$rRNA_fa)
  #write.table(rbind(c(Genome = id, genome_rRNAfeatures)),
  #            file = file.path(outdir, paste0(id, "_rRNA_features.txt")),
  #            row.names = F, col.names = T, sep = "\t", quote = F)
  # 4
  genomesize = genome_genomicfeatures["genomic Total Size"]
  message("Extracting CDS features")
  genome_orffeatures = analyze_orfs(cds_file, genomesize)
  # write.table(rbind(c(Genome = id, genome_orffeatures)),
  #             file = file.path(outdir, paste0(id, "_ORF_features.txt")),
  #             row.names = F, col.names = T, sep = "\t", quote = F)
  # 5
  message("Extracting protein features")
  genome_proteinfeatures = analyze_proteins(proteins_file)
  # write.table(rbind(c(Genome = id, genome_proteinfeatures)),
  #             file = file.path(outdir, paste0(id, "_protein_features.txt")),
  #             row.names = F, col.names = T, sep = "\t", quote = F)
  # 2
  if(tRNA == TRUE) {
    message("Extracting tRNA features")
    trna_results = run.tRNAscan(genome_file = genome_file)
    if(trna_results$found) {
      genome_tRNAfeatures = analyze_tRNA(trna_results$tRNA_lcfa)
      # write.table(rbind(c(Genome = id, genome_tRNAfeatures)),
      #             file = file.path(outdir, paste0(id, "_tRNA_features.txt")),
      #             row.names = F, col.names = T, sep = "\t", quote = F)
    }

    features = data.frame(Genome = id,
                          feature =
                            c(names(genome_genomicfeatures),
                              names(genome_tRNAfeatures),
                              names(genome_orffeatures),
                              names(genome_proteinfeatures)),
                          value =
                            as.numeric(
                              c(genome_genomicfeatures,
                                genome_tRNAfeatures,
                                genome_orffeatures,
                                genome_proteinfeatures))) %>%
      dplyr::as_tibble()
  } else {
    features = data.frame(Genome = id,
                          feature =
                            c(names(genome_genomicfeatures),
                              names(genome_orffeatures),
                              names(genome_proteinfeatures)),
                          value =
                            as.numeric(
                              c(genome_genomicfeatures,
                                genome_orffeatures,
                                genome_proteinfeatures))) %>%
      dplyr::as_tibble()
  }
  return(features)
}

#' @param genome_features genome_features
#' @importFrom dplyr pull as_tibble left_join mutate select filter
#' @returns model.prediction
#' @export run_ogtmodel
run_ogtmodel <- function(genome_features, model = "genomic+tRNA+ORF+protein") {
  message("Running model for optimum T")
  model.select = models %>%
    dplyr::filter(rank == "superkingdom" & clade == "all_species" & model_type == model) %>%
    dplyr::filter(value != 0) %>%
    dplyr::mutate(value = as.numeric(value)) %>%
    dplyr::select(c("feature", "value"))

  model.select.features = model.select %>% dplyr::pull(feature)
  model.select.intercept = model.select %>% dplyr::filter(feature == "intercept") %>% dplyr::pull(value)

  model.sum = genome_features %>%
    dplyr::filter(feature %in% model.select.features) %>%
    dplyr::left_join(model.select, by = c("feature" = "feature")) %>%
    dplyr::mutate(total = `value.x` * `value.y`) %>%
    dplyr::select(total) %>% sum()
  model.prediction = model.sum + model.select.intercept
  return(model.prediction)
}
