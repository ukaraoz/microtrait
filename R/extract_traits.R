#' Apply rules to a set of detected genes
#'
#' @param fasta_file type out_dir
#'
#' @return
#'
#' @export
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
extracttraits <- function(in_file = system.file("extdata/examples/2619619645/in", "2619619645.fna", package = "microtrait", mustWork = TRUE),
                          type = "genomic", save_tempfiles = F, out_dir = getwd()) {
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
    proteins_file = run.prodigal(genome_file = fasta_file, faa_file = tempfile())
    nseq = countseq.fasta(proteins_file)
    tictoc::toc(log = "TRUE")

    tictoc::tic("run.hmmsearch")
    microtrait_domtblout_file = run.hmmsearch(faa_file = proteins_file, hmm = "microtrait")
    dbcan_domtblout_file = run.hmmsearch(faa_file = proteins_file, hmm = "dbcan")
    tictoc::toc(log = "TRUE")
  }
  if(type == "protein") {  # run gene finder
    id = fs::path_file(in_file)
    id = gsub("\\.fna$|\\.faa$|\\.fa$|", "", id, perl = T)

    proteins_file = in_file
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
  result$id = id
  result$norfs = nseq
  result$genes_detected = map.traits.result$genes_detected
  result$genes_detected_table = map.traits.result$genes_detected_table

  result$domains_detected = map.traits.result$domains_detected
  result$rules_asserted = map.traits.result$rules_asserted
  result$all_traits = map.traits.result$all_traits
  result$time_log = tictoc::tic.log()
  saveRDS(result,
          file = file.path(out_dir, paste0(id, ".microtrait.rds")))
  result
}

#' Apply rules to a set of detected genes
#'
#' @param domtblout_file1 domtblout_file2 out_file
#'
#' @return
#'
#' @export
#'
#' @examples
map.traits.fromdomtblout <- function(microtrait_domtblout, dbcan_domtblout) {
  result <- c(call = match.call())

  genes_detected_table = detect.genes.fromdomtblout(microtrait_domtblout)
  genes_detected = dplyr::pull(genes_detected_table, "hmm_name")

  domains_detected = detect.domains.fromdomtblout(dbcan_domtblout)
  rules_asserted = assert.rulesforgenome(genes_detected, domains_detected)

  binary_traits = count.traitsforgenome(rules_asserted, "binary")
  count_traits = count.traitsforgenome(rules_asserted, "count")
  all_traits = dplyr::bind_rows(binary_traits, count_traits)

  result$genes_detected_table = genes_detected_table
  result$genes_detected = genes_detected
  result$domains_detected = domains_detected
  result$rules_asserted = rules_asserted
  result$all_traits = all_traits
  result
}
