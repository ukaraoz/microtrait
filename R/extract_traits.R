#' Apply rules to a set of detected genes
#'
#' @param fasta_file type out_dir
#'
#' @return
#' @import fs dplyr tictoc
#' @importFrom assertthat assert_that
#' @export extract.traits
#'
#' @examples
#' \dontrun{in_file = system.file("extdata/examples/2619619645/in", "2619619645.genes.faa", package = "microtrait", mustWork = TRUE)
#' in_file = list(system.file("extdata/examples/2619619645/out", "2619619645.genes.faa.microtrait.domtblout", package = "microtrait", mustWork = TRUE),
#'                system.file("extdata/examples/2619619645/out", "2619619645.genes.faa.dbcan.domtblout", package = "microtrait", mustWork = TRUE))
#' }
#' \dontrun{
#' in_file = list(system.file("extdata/examples/2619619645/out", "2619619645.genes.faa.microtrait.domtblout", package = "microtrait", mustWork = TRUE),
#'                system.file("extdata/examples/2619619645/out", "2619619645.genes.faa.dbcan.domtblout", package = "microtrait", mustWork = TRUE))
#' }
#' type = "domtblout"
extract.traits <- function(in_file = system.file("extdata/genomic", "2503283023.fna", package = "microtrait", mustWork = TRUE),
                           out_dir = system.file("extdata/genomic", package = "microtrait", mustWork = TRUE),
                           type = "genomic", mode = "single", save_tempfiles = F) {
  result <- c(call = match.call())

  tictoc::tic.clearlog()
  tictoc::tic("extract.traits")

  # Check arguments
  if(type %in% c("genomic", "protein")){
    in_file = tryCatch({
      assertthat::assert_that(file.exists(in_file))
      in_file
    },
    error = function(e) {
      message("Input file ", `in_file`, " doesn't exist.")
      print(e)
    }
    )
  }

  type = tryCatch({
    assertthat::assert_that(type %in% c("genomic", "protein", "domtblout"))
    type
  },
  error = function(e) {
    message("Type has to be either `genomic` or `c`  or `domtblout`.")
    print(e)
  }
  )

  if(type == "genomic") {  # run gene finder
    id = fs::path_file(in_file)
    id = gsub("\\.fna$|\\.faa$|\\.fa$|", "", id, perl = T)

    tictoc::tic("run.prodigal")
    fasta_file = in_file
    proteins_file = run.prodigal(genome_file = fasta_file, faa_file = tempfile(), mode = mode)
    nseq = countseq.fasta(proteins_file)
    tictoc::toc(log = "TRUE")

    if(nseq == 0) {
      map.traits.result = list()
      map.traits.result$genes_detected_table = NULL
      map.traits.result$genes_detected = NULL
      map.traits.result$domains_detected = NULL
      map.traits.result$rules_asserted = NULL
      map.traits.result$trait_counts_atgranularity1 = NULL
      map.traits.result$trait_counts_atgranularity2 = NULL
      map.traits.result$trait_counts_atgranularity3 = NULL
      map.traits.result$id = id
      map.traits.result$norfs = nseq
      map.traits.result$time_log = tictoc::tic.log()
      returnList = list(microtrait_result = map.traits.result, rds_file = NULL)
      return(returnList)
    } else {
      tictoc::tic("run.hmmsearch")
      microtrait_domtblout_file = run.hmmsearch(faa_file = proteins_file, hmm = "microtrait")
      dbcan_domtblout_file = run.hmmsearch(faa_file = proteins_file, hmm = "dbcan")
      tictoc::toc(log = "TRUE")
    }
  }
  if(type == "protein") {  # skip gene finder
    id = fs::path_file(in_file)
    id = gsub("\\.fna$|\\.faa$|\\.fa$|", "", id, perl = T)

    proteins_file = in_file
    nseq = countseq.fasta(proteins_file)
    tictoc::tic("run.hmmsearch")
    microtrait_domtblout_file = run.hmmsearch(faa_file = proteins_file, hmm = "microtrait")
    dbcan_domtblout_file = run.hmmsearch(faa_file = proteins_file, hmm = "dbcan")
    tictoc::toc(log = "TRUE")
  }
  if(type == "domtblout") {  # run gene finder
    id = fs::path_file(in_file[[1]])
    id = gsub("\\.fna$|\\.faa$|\\.fa$|", "", id, perl = T)

    microtrait_domtblout_file = in_file[[1]]
    dbcan_domtblout_file = in_file[[2]]
  }

  tictoc::tic("read.domtblout")
  # microtrait_domtblout_file = "/Users/ukaraoz/Work/microtrait/code/microtrait/inst/extdata/examples/2619619645/out/2619619645.genes.faa.microtrait.domtblout"
  microtrait_domtblout <- read.domtblout(microtrait_domtblout_file)
  # dbcan_domtblout_file = "/Users/ukaraoz/Work/microtrait/code/microtrait/inst/extdata/examples/2619619645/out/2619619645.genes.faa.dbcan.domtblout"
  dbcan_domtblout <- read.domtblout(dbcan_domtblout_file)
  tictoc::toc(log = "TRUE")

  tictoc::tic("map.traits.fromdomtblout")
  map.traits.result = map.traits.fromdomtblout(microtrait_domtblout, dbcan_domtblout)
  tictoc::toc(log = "TRUE")

  tictoc::toc(log = "TRUE")

  if(save_tempfiles == T) {
    file.copy(proteins_file,
              file.path(out_dir, paste0(id, ".prodigal.faa")))
    file.copy(microtrait_domtblout_file,
              file.path(out_dir, paste0(id, ".microtrait.domtblout")))
    file.copy(dbcan_domtblout_file,
              file.path(out_dir, paste0(id, ".dbcan.domtblout")))
  }
  map.traits.result$id = id
  map.traits.result$norfs = nseq
  map.traits.result$time_log = tictoc::tic.log()
  rds_file = file.path(out_dir, paste0(id, ".microtrait.rds"))
  saveRDS(map.traits.result,file = rds_file)
  returnList = list(microtrait_result = map.traits.result, rds_file = rds_file)
  return(returnList)
}

#' Apply rules to a set of detected genes
#'
#' @param domtblout_file1 domtblout_file2 out_file
#' @import dplyr
#' @return
#'
#' @export map.traits.fromdomtblout
map.traits.fromdomtblout <- function(microtrait_domtblout, dbcan_domtblout) {
  result <- c(call = match.call())

  genes_detected_table = detect.genes.fromdomtblout(microtrait_domtblout)
  genes_detected = dplyr::pull(genes_detected_table, "hmm_name")

  domains_detected = detect.domains.fromdomtblout(dbcan_domtblout)
  rules_asserted = assert.rulesforgenome(genes_detected, domains_detected)

  trait_counts_atgranularity1 = count.traitsforgenome(rules_asserted, trait_granularity = 1)
  trait_counts_atgranularity2 = count.traitsforgenome(rules_asserted, trait_granularity = 2)
  trait_counts_atgranularity3 = count.traitsforgenome(rules_asserted, trait_granularity = 3)

  result$genes_detected_table = genes_detected_table
  result$genes_detected = genes_detected
  result$domains_detected = domains_detected
  result$rules_asserted = rules_asserted
  result$trait_counts_atgranularity1 = trait_counts_atgranularity1
  result$trait_counts_atgranularity2 = trait_counts_atgranularity2
  result$trait_counts_atgranularity3 = trait_counts_atgranularity3
  result
}
