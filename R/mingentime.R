#' predictGrowth from codon usage bias
#'
#' @param cds_file proteins_file
#' @importFrom Biostrings readDNAStringSet
#' @importFrom dplyr pull
#' @import gRodon
#' @returns
#' predicted minimum generation time.
#' @export run.predictGrowth
run.predictGrowth <- function(cds_file, proteins_file, mode = "full", temperature = "none", training_set = "vs", depth_of_coverage = NULL, fragments = FALSE) {
  tictoc::tic.clearlog()
  tictoc::tic("predictGrowth")
  message("Running predictGrowth for ", fs::path_file(cds_file), " and ", fs::path_file(proteins_file))

  result <- c(call = match.call())

  ribosomal_domtblout_file = run.hmmsearch(faa_file = proteins_file, hmm = "ribosomal")
  ribosomal_domtblout <- read.domtblout(ribosomal_domtblout_file)

  genes <- Biostrings::readDNAStringSet(cds_file)
  names(genes) = sub(" .*", "", names(genes))
  highlyexpressed_genes = ribosomal_domtblout %>% dplyr::pull(`gene_name`)
  ribosomal_genes = ribosomal_domtblout %>% dplyr::select(c("gene_name", "hmm_name"))
  highlyexpressed_logical = tryCatch(names(genes) %in% (highlyexpressed_genes),error=function(e) NULL)

  if(!is.null(highlyexpressed_logical)) {
    nhighlyexpressed = length(which(highlyexpressed_logical == TRUE))
    temp = myTryCatch(gRodon::predictGrowth(genes,
                                            highlyexpressed_logical,
                                            mode = mode,
                                            temperature = temperature,
                                            training_set = training_set,
                                            depth_of_coverage = NULL,
                                            fragments = FALSE))
    returnList = temp$value
    returnList$warning = temp$warning
    returnList$error = temp$error
    returnList$genes_file = basename(cds_file)
    returnList$highlyexpressed_genes = highlyexpressed_genes
    returnList$nhighlyexpressed = length(highlyexpressed_genes)
    returnList$ribosomal_domtblout_file = ribosomal_domtblout_file
  } else {
    returnList$CUBHE = NA
    returnList$ConsistencyHE = NA
    returnList$CPB = NA
    returnList$FilteredSequences = NA
    returnList$d = NA
    returnList$LowerCI = NA
    returnList$UpperCI = NA
    returnList$genes_file = basename(cds_file)
    returnList$highlyexpressed_genes = highlyexpressed_genes
    returnList$nhighlyexpressed = length(highlyexpressed_genes)
    returnList$ribosomal_domtblout_file = ribosomal_domtblout_file
    returnList$warning = NA
    returnList$error = NA
  }
  tictoc::toc(log = "TRUE")
  returnList
}
