#' Count prevalence of traits, rules, or hmms
#'
#' @param genomes_matrix
#' @param type
#' @import dplyr
#' @return results
#'
#' @export run.parallel
compute.prevalence <- function(feature_matrix, type) {
  if(type == "hmm_matrix") {
    prevalence = feature_matrix %>% tibble::as_tibble() %>%
      dplyr::select(c("id", intersect(colnames(feature_matrix), hmms_fromrules %>% pull(`microtrait_hmm-name`)))) %>%
      tidyr::pivot_longer(!id, names_to = "hmm", values_to = "presence") %>%
      dplyr::select(-id) %>%
      dplyr::count(hmm, presence,.drop=FALSE) %>%
      dplyr::group_by(hmm) %>%
      dplyr::mutate(percent = n/sum(n)*100) %>%
      dplyr::filter(`presence` == 0) %>%
      dplyr::mutate(`percent` = 100-`percent`) %>%
      dplyr::select(c("hmm", "percent")) %>%
      dplyr::inner_join(hmms_fromrules, by = c("hmm" = "microtrait_hmm-name"))
  }

  if(type == "rule_matrix") {
    prevalence = feature_matrix %>% tibble::as_tibble() %>%
      dplyr::select(c("id", rules %>% pull(`microtrait_rule-name`))) %>%
      tidyr::pivot_longer(!id, names_to = "rule", values_to = "presence") %>%
      dplyr::select(-id) %>%
      dplyr::count(rule, presence,.drop=FALSE) %>%
      dplyr::group_by(rule) %>%
      dplyr::mutate(percent = n/sum(n)*100) %>%
      dplyr::filter(`presence` == 0) %>%
      dplyr::mutate(`percent` = 100-`percent`) %>%
      dplyr::select(c("rule", "percent")) %>%
      dplyr::inner_join(rules, by = c("rule" = "microtrait_rule-name"))
  }

  if(type == "trait_matrixatgranularity3") {
    prevalence = feature_matrix %>% tibble::as_tibble() %>%
      dplyr::select(c("id",
                      traits_listbygranularity[[3]] %>%
                        dplyr::filter(`microtrait_trait-type` == "binary") %>%
                        dplyr::pull(`microtrait_trait-name`))) %>%
      dplyr::mutate_at(vars(!starts_with("id")),
                       funs(case_when(. >= 1 ~ 1,TRUE ~ 0))) %>%
      tidyr::pivot_longer(!id, names_to = "trait", values_to = "presence") %>%
      dplyr::select(-id) %>%
      dplyr::count(trait, presence,.drop=FALSE) %>%
      dplyr::group_by(trait) %>%  methaneoxidation_genes = c("mmoX", "mmoY", "mmoZ", "mmoC", "mmoD", "pmoA-amoA", "pmoB-amoB", "pmoC-amoC")

      dplyr::mutate(percent = n/sum(n)*100) %>%
      dplyr::filter(`presence` == 0) %>%
      dplyr::mutate(`percent` = 100-`percent`) %>%
      dplyr::select(c("trait", "percent"))
  }
  prevalence
}
