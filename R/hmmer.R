#' Read a file created through the '----tblout' option
#'
#' @param file (Required). hmmscan output file to be read.
#' @importFrom checkmate expect_file expect_choice expect_number
#' @return data frame of class tibble
#' @export read.tblout
read.tblout <- function(file){
  checkmate::expect_file(file)
  parse.hmmeroutput(file, 'tblout')
}

#' Read a file created through the '----domtblout' option
#'
#' @param file (Required). hmmscan output file to be read.
#' @importFrom checkmate expect_file expect_choice expect_number
#' @return data frame of class tibble
#' @export read.domtblout
read.domtblout <- function(file){
  checkmate::expect_file(file)
  parse.hmmeroutput(file, 'domtblout')
}

#' @import readr
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
parse.hmmeroutput <- function(file, type){
  #checkmate::expect_file(file)
  #checkmate::expect_choice(type,
  #                         choices = c("tblout", "domtblout"),
  #                         info = "hmmer output type is invalid.")
  #file = system.file("extdata/out", "2619619645.genes.faa.dbcan.domtblout", package = "microtrait", mustWork = TRUE)
  col_types <- if(type == 'tblout'){
    readr::cols(
      gene_name         = readr::col_character(),
      gene_accession    = readr::col_character(),
      hmm_name          = readr::col_character(),
      hmm_accession     = readr::col_character(),
      hmm_evalue        = readr::col_double(),
      hmm_score         = readr::col_double(),
      hmm_bias          = readr::col_double(),
      best_gene_evalue  = readr::col_double(),
      best_gene_score   = readr::col_double(),
      best_gene_bis     = readr::col_double(),
      gene_number_exp   = readr::col_double(),
      gene_number_reg   = readr::col_integer(),
      gene_number_clu   = readr::col_integer(),
      gene_number_ov    = readr::col_integer(),
      gene_number_env   = readr::col_integer(),
      gene_number_dom   = readr::col_integer(),
      gene_number_rep   = readr::col_integer(),
      gene_number_inc   = readr::col_character(),
      description       = readr::col_character()
    )
  } else if(type == 'domtblout'){
    readr::cols(
      gene_name         = readr::col_character(),
      gene_accession    = readr::col_character(),
      gene_len          = readr::col_integer(),
      hmm_name          = readr::col_character(),
      hmm_accession     = readr::col_character(),
      hmm_qlen          = readr::col_integer(),
      hmm_evalue        = readr::col_double(),
      hmm_score         = readr::col_double(),
      hmm_bias          = readr::col_double(),
      gene_N            = readr::col_integer(),
      gene_of           = readr::col_integer(),
      gene_cevalue      = readr::col_double(),
      gene_ievalue      = readr::col_double(),
      gene_score        = readr::col_double(),
      gene_bias         = readr::col_double(),
      hmm_from          = readr::col_integer(),
      hmm_to            = readr::col_integer(),
      gene_from         = readr::col_integer(),
      gene_to           = readr::col_integer(),
      env_from          = readr::col_integer(),
      env_to            = readr::col_integer(),
      acc               = readr::col_double()
    )
  }
  N <- length(col_types$cols)

  readr::read_lines(file) %>%
    sub(
      pattern = sprintf("(%s) *(.*)", paste0(rep('\\S+', N-1), collapse=" +")),
      replacement = '\\1\t\\2',
      perl = TRUE
    ) %>%
    paste0(collapse="\n") %>%
    readr::read_tsv(col_names=c('a', 'acc'), comment='#', na='-') %>%
    tidyr::separate(.data$a, head(names(col_types$cols), -1), sep=' +') %>%
    readr::type_convert(col_types=col_types)
}

#' Read a file created through the '----domtblout' option
#'
#' @param file Filename
#' @import dplyr
#' @importFrom checkmate expect_file expect_choice expect_number
#' @return data.frame
#' @export filter.bycoverage.domtblout
filter.bycoverage.domtblout <- function(domtblout, cov_by, cov_threshold = 60){
  checkmate::expect_choice(cov_by,
                           choices = c("gene", "domain"),
                           info = "cov_by type is invalid.")
  checkmate::expect_number(cov_threshold, lower = 0, upper = 100)

  temp = domtblout %>%
    dplyr::group_by(gene_name,hmm_name) %>%
    dplyr::arrange( desc(hmm_score) ) %>%
    # Pick the top 1 value
    dplyr::slice(1) %>%
    dplyr::mutate(cov_gene = ((gene_to-gene_from)/gene_len)*100,
                  cov_domain = ((hmm_to-hmm_from)/hmm_qlen)*100) %>%
    dplyr::select(gene_name, gene_len, hmm_name, gene_score, gene_from, gene_to, cov_gene, hmm_from, hmm_to, cov_domain)
  if(cov_by == "gene"){
    return(temp %>% dplyr::filter(cov_gene >= cov_threshold))
  }
  if(cov_by == "domain"){
    return(temp %>% dplyr::filter(cov_domain >= cov_threshold))
  }
}
