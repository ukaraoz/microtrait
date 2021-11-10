#' Apply rules to a set of detected genes
#'
#' @param domtblout cov_threshold
#' @import dplyr
#' @return
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
  return(hmmhits.gene)
  #return(dplyr::pull(hmmhits.gene, "hmm_name"))
}

#' Apply rules to a set of detected genes
#'
#' @param domtblout cov_threshold
#' @import dplyr
#' @return
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
#' @param rules_asserted trait_granularity
#' @import dplyr
#' @return
#'
#' @export count.traitsforgenome
count.traitsforgenome <- function(rules_asserted, trait_granularity = 3) {
  options( warn = -1 )

  traitcounts = rule2trait %>%
    dplyr::inner_join(rules_asserted, by = c("microtrait_rule-name" = "microtrait_rule-name",
                                             "microtrait_rule-boolean" = "microtrait_rule-boolean")) %>%
    dplyr::filter(`microtrait_rule-asserted` == TRUE)

  # number of dimensions per granularity: 37 (25 + 12), 100 (66+34), 189 (112+77)
  if(trait_granularity == 1) {
    traitcounts = traitcounts %>%
      dplyr::inner_join(traits_listbygranularity[[1]], by = c("microtrait_trait-name1" = "microtrait_trait-name")) %>%
      dplyr::select(c("microtrait_rule-name", "microtrait_rule-boolean", "microtrait_rule-type",
                      "microtrait_rule-substrate", "microtrait_trait-name1", "microtrait_trait-displaynameshort",
                      "microtrait_trait-displaynamelong","microtrait_trait-strategy", "microtrait_trait-granularity",
                      "microtrait_trait-displayorder")) %>%
      dplyr::group_by(`microtrait_trait-name1`, .drop = FALSE) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::rename("microtrait_trait-name" = "microtrait_trait-name1",
                    "microtrait_trait-value" = "n")

    traitcounts = traitcounts %>%
      dplyr::inner_join(traits_listbygranularity[[1]], by = c("microtrait_trait-name" = "microtrait_trait-name")) %>%
      dplyr::mutate(`microtrait_trait-value1` = case_when(`microtrait_trait-type` == "count" ~ `microtrait_trait-value`,
                                                          `microtrait_trait-type` == "binary" & `microtrait_trait-value` >= 1 ~ as.integer(1),
                                                          `microtrait_trait-type` == "binary" & `microtrait_trait-value` == 0 ~ as.integer(0))
      )

  }

  if(trait_granularity == 2) {
    traitcounts = traitcounts %>%
      dplyr::inner_join(traits_listbygranularity[[2]], by = c("microtrait_trait-name2" = "microtrait_trait-name")) %>%
      dplyr::select(c("microtrait_rule-name", "microtrait_rule-boolean", "microtrait_rule-type",
                      "microtrait_rule-substrate", "microtrait_trait-name2", "microtrait_trait-displaynameshort",
                      "microtrait_trait-displaynamelong","microtrait_trait-strategy", "microtrait_trait-granularity",
                      "microtrait_trait-displayorder")) %>%
      dplyr::group_by(`microtrait_trait-name2`, .drop = FALSE) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::rename("microtrait_trait-name" = "microtrait_trait-name2",
                    "microtrait_trait-value" = "n")

    traitcounts = traitcounts %>%
      dplyr::inner_join(traits_listbygranularity[[2]], by = c("microtrait_trait-name" = "microtrait_trait-name")) %>%
      dplyr::mutate(`microtrait_trait-value1` = case_when(`microtrait_trait-type` == "count" ~ `microtrait_trait-value`,
                                                          `microtrait_trait-type` == "binary" & `microtrait_trait-value` >= 1 ~ as.integer(1),
                                                          `microtrait_trait-type` == "binary" & `microtrait_trait-value` == 0 ~ as.integer(0))
      )
  }

  if(trait_granularity == 3) {
    traitcounts = traitcounts %>%
      dplyr::inner_join(traits_listbygranularity[[3]], by = c("microtrait_trait-name3" = "microtrait_trait-name")) %>%
      dplyr::select(c("microtrait_rule-name", "microtrait_rule-boolean", "microtrait_rule-type",
                      "microtrait_rule-substrate", "microtrait_trait-name3", "microtrait_trait-displaynameshort",
                      "microtrait_trait-displaynamelong","microtrait_trait-strategy", "microtrait_trait-granularity",
                      "microtrait_trait-displayorder")) %>%
      dplyr::group_by(`microtrait_trait-name3`, .drop = FALSE) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::rename("microtrait_trait-name" = "microtrait_trait-name3",
                    "microtrait_trait-value" = "n")

    traitcounts = traitcounts %>%
      dplyr::inner_join(traits_listbygranularity[[3]], by = c("microtrait_trait-name" = "microtrait_trait-name")) %>%
      dplyr::mutate(`microtrait_trait-value1` = case_when(`microtrait_trait-type` == "count" ~ `microtrait_trait-value`,
                                                          `microtrait_trait-type` == "binary" & `microtrait_trait-value` >= 1 ~ as.integer(1),
                                                          `microtrait_trait-type` == "binary" & `microtrait_trait-value` == 0 ~ as.integer(0))
      )
  }
  traitcounts
}

#' Apply rules to a set of detected genes
#'
#' @param genes domains cov_threshold
#' @import dplyr
#' @return
assert.rulesforgenome <- function(genes, domains, cov_threshold = 80){
  #genes_detected = detect_genesfromdomtblout(domtblout_file1)
  #domains_detected = detect_domainsfromdomtblout(domtblout_file2)
  #domains_detected = sub(".hmm", "", domains_detected)

  rules_asserted = assert.rules(c(genes, domains))
  return(rules_asserted)
}


#' Apply rules to a set of detected genes
#'
#' @param genes gene list
#' @import dplyr stringr
#' @return
assert.rules <- function(hmms){
  #profile = set.geneprofile(genes)
  #rules_genestatus = unlist(lapply(rules[, "microtrait_rule-booleanunwrapped"],
  #                                 FUN=stringr::str_replace_all,
  #                                 profile))
  #names(rules_genestatus) = dplyr::pull(rules, "microtrait_rule-name")
  #
  #rules_evaluated = unlist(lapply(rules_genestatus,
  #                                function(x) eval(parse(text=paste("as.logical(", x, ")", sep = "")))))
  #rules_evaluated0 = dplyr::bind_cols(`rule-name` = names(rules_evaluated),
  #                                   `rule-boolean` = rules[, "microtrait_rule-boolean"],
  #                                   `rule-genestatus` = rules_genestatus,
  #                                   `rule-asserted` = rules_evaluated)
  profile = set.geneprofile(hmms)
  rules = rules %>%
    dplyr::mutate(`microtrait_rule-booleanunwrapped_set` = stringr::str_replace_all(`microtrait_rule-booleanunwrapped`, profile))
  rules_genestatus = rules %>% pull(`microtrait_rule-booleanunwrapped_set`)
  names(rules_genestatus) = rules %>% pull(`microtrait_rule-name`)
  rules_evaluated = unlist(lapply(rules_genestatus,
                                  function(x) eval(parse(text=paste("as.logical(", x, ")", sep = "")))))
  rules_evaluated = tibble(`microtrait_rule-name` = names(rules_evaluated),
                           `microtrait_rule-asserted` = rules_evaluated) %>%
    dplyr::left_join(rules, by = c("microtrait_rule-name" = "microtrait_rule-name")) %>%
    dplyr::select(c("microtrait_rule-name",
                    "microtrait_rule-boolean",
                    "microtrait_rule-booleanunwrapped_set",
                    "microtrait_rule-asserted"))
  return(rules_evaluated)
}

#' Set gene profiles using detected genes
#'
#' @param genes gene list
#'
#' @return
set.geneprofile <- function(genes){
  microtrait_hmmnames = hmms_fromrules %>% pull(`microtrait_hmm-name`)
  microtrait_hmmnames = add_quotes(microtrait_hmmnames)
  profile = vector(mode = "character", length = length(microtrait_hmmnames))
  profile[1:length(profile)] = "0"
  #names(profile) = paste("'", genes, "'", sep = "")
  names(profile) = microtrait_hmmnames
  profile[add_quotes(genes)] = "1"
  names(profile) = gsub("\\+", "\\\\+", names(profile))
  names(profile) = gsub("\\(", "\\\\(", names(profile))
  names(profile) = gsub("\\)", "\\\\)", names(profile))
  return(profile)
}
