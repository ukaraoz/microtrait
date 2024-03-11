#' Split DNA sequences in a fasta file into single sequence files.
#'
#' @param fasta_file FASTA filename.
#' @param suffix a word to add to each file created. Default is "temp".
#' @param outDir directory to write the output file. Default is tempdir()
#' @param ncores number of cores to use
#'
#' @importFrom Biostrings readBStringSet
#' @importFrom Biostrings writeXStringSet
#' @return a vector of filename names created
#'
#' @export splitSeqsToFiles
splitSeqsToFiles <- function(fasta_file = "/Users/ukaraoz/Work/microtrait/code/inst/extdata/GIG_1C_0-10cm_R1-3.slf.fa",
                             suffix= "temp", outDir = tempdir(), ncores = 4) {
  message("Splitting file ", fasta_file)
  totalSeqs = length(fasta.seqlengths(fasta_file, use.names = FALSE))
  chunks = 1
  starts = seq(0, totalSeqs, by = chunks) ## create chunks of starts

  split.files = parallel::mclapply(1:(length(starts)-1),
                                   function(i) {
                                     query.tmp <- Biostrings::readBStringSet(fasta_file, nrec = 1, skip = starts[i])
                                     #filename.out = file.path(outDir, paste(basename(fasta_file), starts[i], runif(1), suffix, sep = "."))
                                     filename.out = file.path(outDir, paste(names(query.tmp), suffix, "fa", sep = "."))
                                     Biostrings::writeXStringSet(query.tmp, filepath = filename.out, format = "fasta")
                                     filename.out
                                   },
                                   mc.cores = ncores)
  return(unlist(split.files))
}

#' Add metadata
#'
#' @param genomeset_results genomeset_results
#' @param genome_metadata genome_metadata
#' @param genome_metadata_idcol genome_metadata_idcol
#' @importFrom dplyr bind_rows filter left_join mutate rename_with select
#' @return results
#'
#' @export add.metadata
add.metadata <- function(genomeset_results, genome_metadata, genome_metadata_idcol) {
  results = mylist <- sapply(names(genomeset_results),function(x) NULL)
  genome_metadata = genome_metadata %>%
    dplyr::rename_with(~ sub(" ", "_", .x), starts_with("NCBI")) %>%
    dplyr::rename_with(~ sub(" ", "_", .x), starts_with("Ecosystem"))

  for(i in 1:length(genomeset_results)) {
    results[[i]] = genomeset_results[[i]] %>%
      dplyr::left_join(genome_metadata, by = c("id" = genome_metadata_idcol))
      #dplyr::left_join(gp_results, by = c("id" = "id")) %>%
      #dplyr::select(-c("sdgentime", "nHEG", "nNonHEG")) %>%
  }
  results
}

#' Combine microtrait outputs for multiple genomes into trait matrices (genomes x traits)
#'
#' @param rds_files rds_files
#' @param ncores ncores
#' @return genomeset_results
#'
#' @export make.genomeset.results
make.genomeset.results <- function(rds_files, ids = NULL, growthrate = T, optimumT = T, ncores = 1) {
  message("Synthesizing traits at granularity 1")
  trait_matrixatgranularity1 = combine.results(rds_files, type = "trait_atgranularity1", ids = ids, ncores = ncores) %>% tibble::as_tibble()

  message("Synthesizing traits at granularity 2")
  trait_matrixatgranularity2 = combine.results(rds_files, type = "trait_atgranularity2", ids = ids, ncores = ncores) %>% tibble::as_tibble()
  message("Synthesizing traits at granularity 3")
  trait_matrixatgranularity3 = combine.results(rds_files, type = "trait_atgranularity3", ids = ids, ncores = ncores) %>% tibble::as_tibble()

  message("Synthesizing hmms results")
  hmm_matrix = combine.results(rds_files, type = "gene", ids = ids, ncores = ncores) %>% tibble::as_tibble()
  message("Synthesizing rule assertions")
  rule_matrix = combine.results(rds_files, type = "rule", ids = ids, ncores = ncores) %>% tibble::as_tibble()

  if(growthrate == T) {
    mingentimes = combine.results(rds_files, type = "growthrate", ids = ids, ncores = ncores) %>% tibble::as_tibble()
    trait_matrixatgranularity1 = trait_matrixatgranularity1 %>%
      dplyr::left_join(mingentimes, by = c("id" = "id"))
    trait_matrixatgranularity2 = trait_matrixatgranularity2 %>%
      dplyr::left_join(mingentimes, by = c("id" = "id"))
    trait_matrixatgranularity3 = trait_matrixatgranularity3 %>%
      dplyr::left_join(mingentimes, by = c("id" = "id"))
    hmm_matrix = hmm_matrix %>%
      dplyr::left_join(mingentimes, by = c("id" = "id"))
    rule_matrix = rule_matrix %>%
      dplyr::left_join(mingentimes, by = c("id" = "id"))
  }
  if(optimumT == T) {
    optimumTs = combine.results(rds_files, type = "optimumT", ids = ids, ncores = ncores) %>% tibble::as_tibble()
    trait_matrixatgranularity1 = trait_matrixatgranularity1 %>%
      dplyr::left_join(optimumTs, by = c("id" = "id"))
    trait_matrixatgranularity2 = trait_matrixatgranularity2 %>%
      dplyr::left_join(optimumTs, by = c("id" = "id"))
    trait_matrixatgranularity3 = trait_matrixatgranularity3 %>%
      dplyr::left_join(optimumTs, by = c("id" = "id"))
    hmm_matrix = hmm_matrix %>%
      dplyr::left_join(optimumTs, by = c("id" = "id"))
    rule_matrix = rule_matrix %>%
      dplyr::left_join(optimumTs, by = c("id" = "id"))
  }

  genomeset_results = list(trait_matrixatgranularity1 = trait_matrixatgranularity1,
                           trait_matrixatgranularity2 = trait_matrixatgranularity2,
                           trait_matrixatgranularity3 = trait_matrixatgranularity3,
                           hmm_matrix = hmm_matrix,
                           rule_matrix = rule_matrix)
  genomeset_results
}

#' Add metadata
#'
#' @param rds_dir rds_dir
#' @param ncores ncores
#' @return genomeset_results
#'
#' @export make.genomeset.results
write.genomeset.results <- function(genomeset_results, out_dir, datasetid) {
  objects = c("trait_matrixatgranularity1", "trait_matrixatgranularity2", "trait_matrixatgranularity3",
              "hmm_matrix", "rule_matrix")
  for(i in 1:length(objects)) {
    feature_matrix = genomeset_results[[objects[i]]]
    write.table(feature_matrix,
                file = file.path(out_dir, paste0(datasetid, ".", objects[i], ".xls")),
                row.names = F, col.names = T, sep = "\t", quote = F)
  }
}


#' Run microtrait for genomes in parallel using mclapply
#'
#' @param fna_files cache_dir out_dir
#' @import tictoc parallel doParallel foreach
#' @return results
#'
#' @export extract.traits.parallel
extract.traits.parallel <- function(fa_files,
                                    out_dirs = dirname(fa_files),
                                    ncores = floor(parallel::detectCores()*0.7),
                                    save_tempfiles = F, type = "genomic") {
  tictoc::tic(paste0("Running microtrait for ", length(fa_files), " genomes"))
  results = parallel::mclapply(1:length(fa_files),
                      function(i) {
                         returnList = extract.traits(fa_files[i], out_dirs[i], save_tempfiles = save_tempfiles, type = type)
                         returnList
                       },
                      mc.cores = ncores)

  # alternative way to parallelize, more error-prone in cluster environment due to memory allocation errors
  #cl <- parallel::makeCluster(ncores)
  #doParallel::registerDoParallel(cl)
  #results = foreach(i = 1:length(fa_files), .packages = c("microtrait")) %dopar% {
  #  extract.traits(fa_files[i], save_tempfiles = T, out_dir = out_dir)
  #}
  #parallel::stopCluster(cl)
  tictoc::toc(log = "TRUE")
  results
}

#' Combine microtrait results for multiple genomes
#'
#' @import parallel doParallel tidyr
#'
#' @description
#' `combine.results()` builds trait, rule, and hmm matrices by combining microtrait results
#' for multiple genomes.
#'
#'
#' @param rds_files Full paths to rds files that holds microtrait results.
#' @param ids ids to be assigned to the genomes, default is the rds file name
#' @param type microtrait result to be combined
#' @param ncores Number of cores for computation
#'
#' @import tibble tidyr
#' @return trait, rule, or hmm matrix
combine.results <- function(rds_files, ids = NULL, type, ncores = 1) {
  # reminder for "foreach" way to parallelize
  # rds_files = unlist(purrr::map_depth(results, 1, "rds_file"))
  if (is.null(ids)) {
    ids = sub("\\..*", "", sub(".rds", "", basename(rds_files)), perl = T)
  }
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  results_rowbinded = foreach(i = 1:length(rds_files), .packages = c("dplyr", "tibble"), .export = "fetch.results", .combine = "rbind") %dopar% {
    fetch.results(rds_files[i], ids[i], type = type)
  }
  results_rowbinded = results_rowbinded %>% tibble::as_tibble()
  results_matrix = switch(type,
    trait_atgranularity1 = {results_rowbinded %>% tidyr::spread(`microtrait_trait-name`, `microtrait_trait-value`)},
    trait_atgranularity2 = {results_rowbinded %>% tidyr::spread(`microtrait_trait-name`, `microtrait_trait-value`)},
    trait_atgranularity3 = {results_rowbinded %>% tidyr::spread(`microtrait_trait-name`, `microtrait_trait-value`)},
    # distinct() needed here, multiple hits/genome
    # drop=FALSE to add columns for undetected hmms
    gene = {results_rowbinded %>% tibble::add_column(n = factor(1, levels = c(0,1), ordered = T)) %>% distinct() %>% tidyr::spread(hmm, n, drop = FALSE, fill = 0)},
    rule = {results_rowbinded %>% tidyr::spread(`microtrait_rule-name`, `microtrait_rule-asserted`)},
    mingentime = {results_rowbinded %>% tidyr::spread(`microtrait_trait-name`, `microtrait_trait-value`)},
    optimumT = {results_rowbinded %>% tidyr::spread(`microtrait_trait-name`, `microtrait_trait-value`)}
  )
  parallel::stopCluster(cl)
  # todo: NULLify rownames earlier
  row.names(results_matrix) = NULL
  results_matrix
}

#' Fetch results for a single genome
#'
#' @param rds_file id type
#' @import tibble
#' @returns
#' fetched results from a single genome.
#' @export
fetch.results <- function(rds_file, id, type) {
  temp = readRDS(rds_file)
  if(type == "trait_atgranularity1") {
    result = temp$trait_counts_atgranularity1 %>%
      select(c("microtrait_trait-name", "microtrait_trait-value"))
  }
  if(type == "trait_atgranularity2") {
    result = temp$trait_counts_atgranularity2 %>%
      select(c("microtrait_trait-name", "microtrait_trait-value"))
  }
  if(type == "trait_atgranularity3") {
    result = temp$trait_counts_atgranularity3 %>%
      select(c("microtrait_trait-name", "microtrait_trait-value"))
  }
  if(type == "gene") {
    result = dplyr::bind_rows(temp$genes_detected %>% tibble::as_tibble(),
                              temp$domains_detected %>% tibble::as_tibble()) %>%
             dplyr::rename(hmm = value) %>%
             dplyr::mutate(`hmm` = readr::parse_factor(`hmm`,
                                                       ordered = T,
                                                       levels = hmms_fromrules %>% dplyr::pull(`microtrait_hmm-name`) %>% levels)
                          )
  }
  if(type == "rule") {
    result = temp$rules_asserted %>%
      select(c("microtrait_rule-name", "microtrait_rule-asserted"))
  }
  # the returned variable is name growthrate but it is not a rate, it is mingentime in hours
  if(type == "growthrate") {
    result = data.frame(`microtrait_trait-name` = "mingentime",
                        `microtrait_trait-value` = temp$growthrate_d,
                        check.names = F) %>% tibble::as_tibble()
  }
  if(type == "optimumT") {
    result = data.frame(`microtrait_trait-name` = "optimumT",
                        `microtrait_trait-value` = temp$ogt,
                        check.names = F) %>% tibble::as_tibble()
  }
  result = result %>%
    tibble::add_column(id = as.character(id)) %>%
    dplyr::select(id, everything()) %>% as.data.frame
  return(result)
}

combine.gp.results <- function(files, ids, ncores = 1) {
  parsed = mclapply(files, gp.resultsummary, mc.cores = ncores)
  temp.mingentime <- mclapply(parsed, "[[", 1, mc.cores = ncores)
  temp.sdgentime <- mclapply(parsed, "[[", 2, mc.cores = ncores)
  temp.nHEG <- mclapply(parsed, "[[", 3, mc.cores = ncores)
  temp.nNonHEG <- mclapply(parsed, "[[", 4, mc.cores = ncores)
  temp.OGT <- mclapply(parsed, "[[", 5, mc.cores = ncores)

  results = data.frame(id = ids,
                       mingentime = as.numeric(unlist(temp.mingentime)),
                       sdgentime = as.numeric(unlist(temp.sdgentime)),
                       nHEG = as.numeric(unlist(temp.nHEG)),
                       nNonHEG = as.numeric(unlist(temp.nNonHEG)),
                       OGT = as.numeric(unlist(temp.OGT)),
                       stringsAsFactors = F) %>% tibble::as_tibble()
  results
}

#' Fetch growthpred results for a single genome
#'
#' @param file file
#' @returns
#' fetches growthpred results.
gp.resultsummary <- function(file) {
  temp = readLines(file, n = -1)
  mingentime = gsub("Predicted minimum generation time:  (.*) hours \\+/- (.*)  ", "\\1", temp[1], perl=TRUE)
  sdgentime = gsub("Predicted minimum generation time:  (.*) hours \\+/- (.*)  ", "\\2", temp[1], perl=TRUE)
  nHEG = gsub("Input dataset: (.*) highly expressed genes \\+ (.*) non-highly expressed genes ", "\\1", temp[3], perl=TRUE)
  nNonHEG = gsub("Input dataset: (.*) highly expressed genes \\+ (.*) non-highly expressed genes ", "\\2", temp[3], perl=TRUE)
  OGT = gsub("Input optimal growth temperature \\(Celsius\\): (.*) ", "\\1", temp[4], perl=TRUE)
  result = list(mingentime = mingentime,
                sdgentime = sdgentime,
                nHEG = nHEG,
                nNonHEG = nNonHEG,
                OGT = OGT
  )
  result
}

# files = list.files("/Users/ukaraoz/Work/microtrait/code/inst/extdata/grodon", pattern = "_growth.rds$", full.names = T)
combine.grodon.results <- function(files, ids = gsub("\\..*", "", basename(files), perl = T), ncores = 1) {
  parsed = mclapply(files, readRDS, mc.cores = ncores)
  # mean of the codon usage bias of each highly expressed gene relative to all other genes
  temp.CUBHE <- mclapply(parsed, "[[", 1, mc.cores = ncores)

  #mean of the codon usage bias of each highly expressed gene relative to all other highly expressed genes
  temp.ConsistencyHE <- mclapply(parsed, "[[", 2, mc.cores = ncores)

  #genome-wide codon pair bias
  temp.codonpairbias <- mclapply(parsed, "[[", 3, mc.cores = ncores)

  #number of sequences filtered due to length
  temp.FilteredSequences <- mclapply(parsed, "[[", 4, mc.cores = ncores)

  #estimated minimal doubling time in hours
  temp.mingentime <- mclapply(parsed, "[[", 5, mc.cores = ncores)

  #95% confidence intervals for `d`
  temp.LowerCI <- mclapply(parsed, "[[", 6, mc.cores = ncores)
  temp.UpperCI <- mclapply(parsed, "[[", 7, mc.cores = ncores)
  temp.genes_file <- mclapply(parsed, "[[", 8, mc.cores = ncores)
  temp.highlyexpressed_file <- mclapply(parsed, "[[", 9, mc.cores = ncores)
  temp.nhighlyexpressed <- mclapply(parsed, "[[", 10, mc.cores = ncores)

  results = data.frame(id = ids,
                       CUBHE = as.numeric(unlist(temp.CUBHE)),
                       ConsistencyHE = as.numeric(unlist(temp.ConsistencyHE)),
                       codonpairbias = as.numeric(unlist(temp.codonpairbias)),
                       FilteredSequences = as.numeric(unlist(temp.FilteredSequences)),
                       mingentime = as.numeric(unlist(temp.mingentime)),
                       LowerCI = as.numeric(unlist(temp.LowerCI)),
                       UpperCI = as.numeric(unlist(temp.UpperCI)),
                       genes_file = unlist(temp.genes_file),
                       highlyexpressed_file = unlist(temp.highlyexpressed_file),
                       nhighlyexpressed = as.numeric(unlist(temp.nhighlyexpressed)),
                       stringsAsFactors = F) %>% tibble::as_tibble()
  results
}
