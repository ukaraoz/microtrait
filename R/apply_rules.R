#' Apply rules to a set of detected genes
#'
#' @param genes_detected
#'
#' @return
#'
#' @examples
detect.genes.fromdomtblout <- function(domtblout, cov_threshold = 80){
  #hmmhits = read_domtblout(domtblout_file)
  # just to make sure
  hmmhits.gene = domtblout %>% dplyr::filter(!grepl("^GH", hmm_name))
  hmmhits.gene =
    hmmhits.gene %>%
    filter.bycoverage.domtblout(cov_by = "gene", cov_threshold) %>%
    dplyr::group_by(gene_name, hmm_name) %>%
    dplyr::top_n(1, cov_gene) %>%
    dplyr::group_by(hmm_name) %>%
    dplyr::top_n(1, gene_score)
  return(dplyr::pull(hmmhits.gene, "hmm_name"))
}

#' Apply rules to a set of detected genes
#'
#' @param genes_detected
#'
#' @return
#'
#' @examples
detect.domains.fromdomtblout <- function(domtblout, cov_threshold = 35){
  hmmhits = domtblout %>%
    dplyr::mutate(hmm_name = stringr::str_replace(hmm_name, ".hmm", ""))
  hmmhits.domain = hmmhits
  hmmhits.domain =
    hmmhits.domain %>%
    filter.bycoverage.domtblout(cov_by = "domain", cov_threshold) %>%
    dplyr::group_by(gene_name, hmm_name) %>%
    dplyr::top_n(1, cov_domain)
  return(dplyr::pull(hmmhits.domain, "hmm_name"))
}

#' Apply rules to a set of detected genes
#'
#' @param genes_detected
#'
#' @return
#'
#' @examples
count.traitsforgenome <- function(rules_asserted, trait) {
  options( warn = -1 )
  if(trait == "binary") {
    traitcounts = rule2trait %>%
      dplyr::filter(`trait-type`=="binary") %>%
      dplyr::inner_join(rules_asserted, by = c("rule-name" = "rule-name")) %>%
      dplyr::arrange(`trait-displayorder`) %>%
      dplyr::select(c("trait-display-short", "rule-asserted", "trait-type")) %>%
      dplyr::mutate(`rule-asserted` = as.integer(`rule-asserted`)) %>%
      dplyr::rename(trait=`rule-asserted`, n = `rule-asserted`)
  }

  if(trait == "count") {
    pattern = "count: \"Stress Tolerance: Biofilm Formation: |\"|count: \"Resource Acquisition: Transport: |\"|count: \"Resource Acquisition: Extracellular: |\"|count: \"Stress Tolerance: Osmoprotection: |\""
    traitcounts = rule2trait %>%
      dplyr::filter(grepl("^count: ", `trait-type`)) %>%
      dplyr::mutate(`trait-type` = gsub(pattern, "", `trait-type`)) %>%
      dplyr::inner_join(rules_asserted, by = c("rule-name" = "rule-name")) %>%
      dplyr::filter(`rule-asserted` == TRUE) %>%
      tidyr::separate_rows(`trait-type`, sep = ";") %>%   # at this point, the trait-type is a countable object, i.e. substrate, they should be all in "substrates"
      dplyr::left_join(substrates, by = c("trait-type" = "substrate-name")) %>%
      dplyr::left_join(substrateclasses, by = c("substrate-class" = "substrate-class")) %>%
      dplyr::select(c(`substrate-class`)) %>%
      dplyr::group_by(`substrate-class`, .drop = FALSE) %>%   # key not to drop levels with zero counts
      dplyr::count(.drop = FALSE) %>%   # have the counts, polish
      dplyr::left_join(substrateclasses, by = c("substrate-class" = "substrate-class")) %>%
      dplyr::arrange(`trait-displayorder`) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("trait-display-short", "n", "trait-type"))
  }
  return(traitcounts)
}

#' Apply rules to a set of detected genes
#'
#' @param genes_detected
#'
#' @return
#'
#' @examples
assert.rulesforgenome <- function(genes, domains, cov_threshold = 80){
  #genes_detected = detect_genesfromdomtblout(domtblout_file1)
  #domains_detected = detect_domainsfromdomtblout(domtblout_file2)
  #domains_detected = sub(".hmm", "", domains_detected)

  rules_asserted = assert.rules(c(genes, domains))
  return(rules_asserted)
}


#' Apply rules to a set of detected genes
#'
#' @param genes
#'
#' @return
#'
#' @examples
assert.rules <- function(genes){
  profile = set.geneprofile(genes)
  rules_genestatus = unlist(lapply(rules[, "rule-boolean-unwrapped"],
                                   FUN=stringr::str_replace_all,
                                   profile))
  names(rules_genestatus) = dplyr::pull(rules, "rule-name")

  rules_evaluated = unlist(lapply(rules_genestatus,
                                  function(x) eval(parse(text=paste("as.logical(", x, ")", sep = ""))))
  )
  #rules_evaluated = dplyr::bind_cols(list(names(rules_evaluated), rules_evaluated))

  rules_evaluated = dplyr::bind_cols(`rule-name` = names(rules_evaluated),
                                     `rule-boolean` = rules[, "rule-boolean"],
                                     `rule-genestatus` = rules_genestatus,
                                     `rule-asserted` = rules_evaluated)
  return(rules_evaluated)
}

#' Set gene profiles using detected genes
#'
#' @param genes
#'
#' @return
#'
#' @examples
set.geneprofile <- function(genes){
  hmmnames = add_quotes(microtraithmm_names)
  profile = vector(mode = "character", length = length(hmmnames))
  profile[1:length(profile)] = "0"
  #names(profile) = paste("'", genes, "'", sep = "")
  names(profile) = hmmnames
  profile[add_quotes(genes)] = "1"
  names(profile) = gsub("\\+", "\\\\+", names(profile))
  names(profile) = gsub("\\(", "\\\\(", names(profile))
  names(profile) = gsub("\\)", "\\\\)", names(profile))
  return(profile)
}
