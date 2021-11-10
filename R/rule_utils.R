#' Return genes that are part of a rule.
#'
#' @param rule (Required). microtrait rule.
#'
#' @return rule parts
#'
#' @export rule2parts
rule2parts <- function(rule_boolean) {
  parts = unlist(strsplit(gsub("\\(+|\\)+|\\!", "", rule_boolean, perl = T), " [\\|\\&]+ ", perl = T))
  return(parts)
}

#' Unwraps complex rules.
#'
#' @param rule (Required). microtrait rule.
#' @import dplyr stringr tidyr
#' @return rule parts
#' @export unwrap.rules
unwrap.rules <- function(rule_table) {
  rule_table = rule_table %>%
    dplyr::mutate(`microtrait_rule-namequoted` = paste0("'", `microtrait_rule-name`, "'"))
  rules_parts = unique(unlist(lapply(rule_table %>% pull(`microtrait_rule-boolean`), FUN=rule2parts)))

  # building a lookup table for parts
  lookuptabletemp = data.frame(`microtrait_rule-part` = rules_parts,
                               `microtrait_rule-partlookup` = rules_parts,
                               check.names = F) %>% tidyr::as_tibble()
  # lookup = parts
  # names(lookup) = lookup
  lookuptabletemp = lookuptabletemp %>%
    dplyr::left_join(rule_table, by = c("microtrait_rule-part" = "microtrait_rule-namequoted")) %>%
    dplyr::filter(!is.na(`microtrait_rule-name`)) %>%
    dplyr::mutate(`microtrait_rule-partlookup` = `microtrait_rule-boolean`) %>%
    dplyr::select(c("microtrait_rule-part", "microtrait_rule-partlookup")) %>%
    dplyr::mutate(`microtrait_rule-part` = gsub("\\+", "\\\\+", `microtrait_rule-part`))
  lookuptable = lookuptabletemp %>% dplyr::pull(`microtrait_rule-partlookup`)
  names(lookuptable) = lookuptabletemp %>% dplyr::pull(`microtrait_rule-part`)

  ## parts2lookup = intersect(lookup, rules %>% pull(`rule-name-quoted`))
  ## lookup[parts2lookup] =
  ##   rules %>%
  ##   filter(`rule-name-quoted` %in% parts2lookup) %>%
  ##   select(`rule-boolean`) %>% pull
  ## names(lookup) = gsub("\\+", "\\\\+", names(lookup))

  rule_table = rule_table %>%
    dplyr::mutate(`microtrait_rule-booleanunwrapped` = rule_table %>%
                                                        pull(`microtrait_rule-boolean`) %>%
                                                        stringr::str_replace_all(lookuptable) %>%
                                                        stringr::str_replace_all(lookuptable) %>%
                                                        stringr::str_replace_all(lookuptable) %>%
                                                        stringr::str_replace_all(lookuptable) %>%
                                                        stringr::str_replace_all(lookuptable) %>%
                                                        stringr::str_replace_all(lookuptable)
                 )
  #rules = rules.df %>% tbl_df
  return(rule_table)
}

#' Get genes from rule
#'
#' @param rule (Required). microtrait rule.
#' @import dplyr
#' @return
rulename2ruleboolean <- function(rule) {
  return(unwrap.rules(rules %>% dplyr::filter(`rule_name` %in% rule)) %>%
           dplyr::select(`rule_boolean`) %>% dplyr::pull()
  )
}

#' Get genes from rule
#'
#' @param rule (Required). microtrait rule.
#'
#' @return
#' @export rule2genes
rule2genes <- function(rule) {
  return(removequotes(unique(rule2parts(rulename2ruleboolean(rule)))))
}

removequotes <- function(x){
  return(gsub("'", "", x))
}

add_quotes <- function(string) {
  return(paste("'", string, "'", sep = ""))
}

#' Evaluate a boolean rule.
#'
#' @param boolean
#'
#' @return
evaluate.rule <- function(boolean){
  eval(parse(text=paste("as.logical(", boolean, ")", sep = "")))
}
