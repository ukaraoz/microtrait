## code to prepare `rule2substrate` dataset goes here
library(tidyverse)
library(usethis)

source("./R/rule_utils.R")

microtraithmm_names = readr::read_delim("./data-raw/microtraithmm-names.txt",
                                       delim = "\t", col_names = F
                                      ) %>% pull()
message(paste("Read", length(microtraithmm_names), "microtraithmm names."))

rules = readr::read_delim("./data-raw/rules.txt",
                          delim = "\t",
                          col_types = readr::cols(
                            `rule-name` = readr::col_character(),
                            `rule-boolean` = readr::col_character()
                            )
                          )
message(paste("Read", nrow(rules), "rules."))
rules = unwrap.rules(rules)
hmm_names_fromrules = unique(unlist(lapply(rules %>% pull(`rule-boolean-unwrapped`), FUN=rule2parts)))

#write.table(rules, file = "./data-raw/rules-unwrapped.txt", row.names = F, col.names = T, sep = "\t")
hmm_names_fromrules = removequotes(hmm_names_fromrules)
hmm_names_fromrules_missing = setdiff(hmm_names_fromrules, microtraithmm_names)
message(paste(length(hmm_names_fromrules_missing), " hmm-names after unwrapping the rules don't match microtraithmm_names.\n This should be 0."))
write.table(rules, file = "./data-raw/rules-unwrapped.xls", row.names = F, col.names = T, sep = "\t")
usethis::use_data(rules, microtraithmm_names, overwrite = TRUE, compress = 'xz')

rule2trait = read_delim("./data-raw/rule2trait.txt", delim = "\t", col_names = T,
                            col_types =
                              cols(`rule-name` = col_character(),
                                   `trait-name` = col_character(),
                                   `trait-display-long` = col_character(),
                                   `trait-display-short` = col_character(),
                                   `trait-type` = col_character(),
                                   `trait-category` = col_character(),
                                   `trait-version` = col_character(),
                                   `trait-displayorder` = col_integer()
                                  )
                          )
message(paste("Read", nrow(rule2trait), "rule-to-trait association."))

# make sure all rules are defined
rule2trait_rulesnotdefined =
  setdiff(rule2trait %>% pull(`rule-name`),
          rules %>% pull(`rule-name`))
message(paste(length(rule2trait_rulesnotdefined), " rules not defined.\n This should be 0."))
usethis::use_data(rule2trait, overwrite = TRUE)

# read substrates and substrate classes
substrates = read_delim("data-raw/substrates.txt", delim = "\t", col_names = T,
                        col_types =
                          cols(`substrate-name` = col_factor(),
                               `substrate-subclass1` = col_factor(),
                               `substrate-subclass2` = col_factor(),
                               `substrate-subclass3` = col_factor(),
                               `substrate-class` = col_factor(),
                               `substrate-version` = col_factor()
                          )) %>%
                        dplyr::filter(`substrate-version` == "production") %>%
                        droplevels
substrateclasses = read_delim("data-raw/substrate-class.txt", delim = "\t", col_names = T,
                              col_types =
                                cols(`substrate-class` = col_factor(),
                                     `trait-display-short` = col_factor(),
                                     `trait-type` = col_factor(),
                                     `trait-version` = col_factor(),
                                     `trait-displayorder` = col_integer()
                                )) %>%
                  dplyr::filter(`trait-version` == "production") %>%
                  droplevels
usethis::use_data(substrates, substrateclasses, overwrite = TRUE)


###
#traits_transport = rules %>%
#  dplyr::inner_join(rule2substrate, by = c("rule_name" = "rule_name")) %>%
#  dplyr::filter(`trait_type`=="transport") %>%
#  tidyr::separate_rows(`substrate_name`, sep = ";") %>%
#  dplyr::inner_join(substrates, by = c("substrate_name" = "substrate_name")) %>%
#  dplyr::select(c(`substrate_class`)) %>%
#  dplyr::rename(trait=substrate_class) %>%
#  unique()
#
#traits_nontransport = rules %>%
#  dplyr::inner_join(rule2substrate, by = c("rule_name" = "rule_name")) %>%
#  dplyr::filter((`rule_class`=="adhesion" |
#                     `rule_class` == "biofilm" |
#                     `rule_class` == "chemotaxis:flagella" |
#                     `rule_class` == "chemotaxis:quorum sensing" |
#                     `rule_class` == "chemotaxis:aerotaxis" |
#                     `rule_class` == "chemotaxis:nutrient sensing" |
#                     `rule_class` == "osmoprotection" |
#                     `rule_class` == "phosphotransferase" |
#                     `rule_class`=="motility" |
#                     `rule_class`=="secretion")) %>%
#  dplyr::select(c(`rule_class`)) %>%
#  dplyr::filter(!(`rule_class` %in% c("uptake", "complex carbohydrate degradation", "drug/toxin resistance", "efflux/export", "metabolic"))) %>%
#  dplyr::rename(trait=rule_class) %>%
#  unique()
#
#traits_GH = rules %>%
#  dplyr::inner_join(rule2substrate, by = c("rule_name" = "rule_name")) %>%
#  #dplyr::filter(`rule_boolean.x` ==  "TRUE" & `rule_class`=="complex carbohydrate degradation") %>%
#  dplyr::filter(`rule_class`=="complex carbohydrate degradation") %>%
#  tidyr::separate_rows(`substrate_name`, sep = ";") %>%
#  dplyr::inner_join(substrates, by = c("substrate_name" = "substrate_name")) %>%
#  dplyr::select(`substrate_name`) %>%
#  dplyr::rename(trait=substrate_name) %>%
#  unique()
#
#traits_metabolic = rules %>%
#  dplyr::inner_join(rule2substrate, by = c("rule_name" = "rule_name")) %>%
#  dplyr::filter(`trait_type`=="metabolic") %>%
#  dplyr::inner_join(trait2displaymetabolic, by = c("rule_name" = "rule_name")) %>%
#  dplyr::select(`trait_displayname`) %>%
#  dplyr::rename(trait=trait_displayname) %>%
#  unique()
#
#traits_transport = traits_transport%>%as.data.frame%>%pull(1)%>%as.character()
#traits_nontransport = traits_nontransport%>%as.data.frame%>%pull(1)%>%as.character()
#traits_GH = traits_GH%>%as.data.frame%>%pull(1)%>%as.character()
#traits_metabolic = traits_metabolic%>%as.data.frame%>%pull(1)%>%as.character()
#alltraits = c(traits_transport,
#                  traits_nontransport,
#                  traits_GH,
#                  traits_metabolic)
#usethis::use_data(traits_metabolic,
#                  traits_transport,
#                  traits_nontransport,
#                  traits_GH,
#                  alltraits, overwrite = TRUE)

