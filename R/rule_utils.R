#' Return genes that are part of a rule.
#'
#' @param rule (Required). microtrait rule.
#'
#' @return rule parts
#'
#' @export
#'
#' @examples
rule2parts <- function(rule_boolean) {
  parts = unlist(strsplit(gsub("\\(+|\\)+", "", rule_boolean, perl = T), " [\\|\\&]+ ", perl = T))
  return(parts)
}

#'
#'
#' @param rule (Required). microtrait rule.
#'
#' @return rule parts
#'
#' @examples
unwrap.rules <- function(rules) {
  rules = rules %>%
    dplyr::mutate(`rule-name-quoted` = paste("'", rules%>%pull('rule-name'), "'", sep = ""))
  parts = unique(unlist(lapply(rules %>% pull(`rule-boolean`), FUN=rule2parts)))

  lookup = parts
  names(lookup) = lookup
  parts2lookup = intersect(lookup, rules %>% pull(`rule-name-quoted`))
  lookup[parts2lookup] =
    rules %>%
    filter(`rule-name-quoted` %in% parts2lookup) %>%
    select(`rule-boolean`) %>% pull
  names(lookup) = gsub("\\+", "\\\\+", names(lookup))

  rules = rules %>%
    dplyr::mutate(`rule-boolean-unwrapped` =
                    rules %>% pull(`rule-boolean`) %>%
                    stringr::str_replace_all(lookup) %>%
                    stringr::str_replace_all(lookup) %>%
                    stringr::str_replace_all(lookup) %>%
                    stringr::str_replace_all(lookup)
    )
  #rules = rules.df %>% tbl_df
  return(rules)
}

#' Get genes from rule
#'
#' @param
#'
#' @return
#'
#' @examples
rulename2ruleboolean <- function(rule) {
  return(unwrap.rules(rules %>% dplyr::filter(`rule_name` %in% rule)) %>%
           dplyr::select(`rule_boolean`) %>% dplyr::pull()
  )
}

#' Get genes from rule
#'
#' @param
#'
#' @return
#' @export
#'
#' @examples
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
#'
#' @examples
evaluate.rule <- function(boolean){
  eval(parse(text=paste("as.logical(", boolean, ")", sep = "")))
}
