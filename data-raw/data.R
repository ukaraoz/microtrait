## code to prepare `data` dataset goes here
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(usethis)

source("./R/rule_utils.R")

# 1. Read information about microtrait hmms
hmms = readr::read_delim("./data-raw/microtrait_hmm.txt",
                                       delim = "\t",
                                       col_types = readr::cols(`microtrait_hmm-name` = readr::col_character(),
                                                               `microtrait_hmm-dbxref_kegg` = readr::col_character(),
                                                               `microtrait_hmm-dbxref_ec` = readr::col_character(),
                                                               `microtrait_hmm-dbxref_TC` = readr::col_character(),
                                                               `microtrait_hmm-description` = readr::col_character(),
                                                               `microtrait_hmm-notes` = readr::col_character()
                                                              )
                                      )
hmms_cazy = hmms %>% dplyr::filter(stringr::str_detect(`microtrait_hmm-name`, "^GH") == TRUE)
message(paste("Read", nrow(hmms), "microtraithmms."))

# 2. Read information about microtrait rules
rules = readr::read_delim("./data-raw/microtrait_rules.txt",
                                         delim = "\t",
                                         col_types = readr::cols(`microtrait_rule-name` = readr::col_factor(ordered = FALSE),
                                                                 `microtrait_rule-boolean` = readr::col_character()
                                                                )
                                        )
message(paste("Read", nrow(rules), "rules."))
# unwrap rules
rules = unwrap.rules(rules)
# find which hmms appear in the rules
hmm_names_fromrules = unique(unlist(lapply(rules %>% pull(`microtrait_rule-booleanunwrapped`), FUN=rule2parts)))
hmm_names_fromrules = removequotes(hmm_names_fromrules)
hmm_names_fromrules = hmm_names_fromrules %>% tibble::as_tibble() %>%
  dplyr::rename(`microtrait_hmm-name` = value)
# This should be empty if properly unwrapped
#hmm_names_fromrules %>%
#  anti_join(hmms, by = c("microtrait_hmm-name" = "microtrait_hmm-name"))
hmm_names_fromrules = hmm_names_fromrules %>%
  dplyr::left_join(hmms, by = c("microtrait_hmm-name" = "microtrait_hmm-name"))

# 3. read traits
traits = readr::read_delim("./data-raw/microtrait_traits.txt",
                           delim = "\t", col_names = T,
                           col_types = readr::cols(`microtrait_trait-name` = readr::col_character(),
                                                   `microtrait_trait-displaynameshort` = readr::col_character(),
                                                   `microtrait_trait-displaynamelong` = readr::col_character(),
                                                   `microtrait_trait-strategy` = readr::col_character(),
                                                   `microtrait_trait-type` = readr::col_character(),
                                                   `microtrait_trait-granularity` = readr::col_factor(levels = c("1","2","3"), ordered = TRUE),
                                                   `microtrait_trait-version` = readr::col_character(),
                                                   `microtrait_trait-displayorder` = readr::col_integer()
                                                  )
                          ) %>%
        dplyr::arrange(`microtrait_trait-strategy`, `microtrait_trait-granularity`, `microtrait_trait-displayorder`)

## split by granularity
traits_granularity_levels = c(1,2,3)
traits_listbygranularity = vector("list", length(traits_granularity_levels))
for(i in 1:length(traits_granularity_levels)) {
  traits_atgranularity = traits %>%
    dplyr::filter(`microtrait_trait-granularity` == traits_granularity_levels[i]) %>%
    dplyr::arrange(`microtrait_trait-displayorder`)
  traits_atgranularity = traits_atgranularity %>%
    dplyr::mutate(`microtrait_trait-name` = readr::parse_factor(`microtrait_trait-name`,
                                                                ordered = T,
                                                                levels = traits_atgranularity %>%
                                                                  dplyr::arrange(`microtrait_trait-displayorder`) %>%
                                                                  dplyr::pull(`microtrait_trait-name`)),
                  `microtrait_trait-displaynameshort` = readr::parse_factor(`microtrait_trait-displaynameshort`,
                                                                            ordered = T,
                                                                            levels = traits_atgranularity %>%
                                                                              dplyr::arrange(`microtrait_trait-displayorder`) %>%
                                                                              dplyr::pull(`microtrait_trait-displaynameshort`)),
                  `microtrait_trait-displaynamelong` = readr::parse_factor(`microtrait_trait-displaynamelong`,
                                                                            ordered = T,
                                                                            levels = traits_atgranularity %>%
                                                                              dplyr::arrange(`microtrait_trait-displayorder`) %>%
                                                                              dplyr::pull(`microtrait_trait-displaynamelong`))
                )
  traits_listbygranularity[[i]] = traits_atgranularity
}

# 4. read substrate to rule
substrate2rule = readr::read_delim("./data-raw/microtrait_substrate2rule.txt",
                                   delim = "\t", col_names = T,
                                   col_types = readr::cols(`microtrait_substrate-name` = readr::col_character(),
                                                           `microtrait_substrate-subclass1` = readr::col_character(),
                                                           `microtrait_substrate-subclass2` = readr::col_character(),
                                                           `microtrait_substrate-subclass3` = readr::col_character(),
                                                           `microtrait_trait-name1` = readr::col_character(),
                                                           `microtrait_trait-name2` = readr::col_character(),
                                                           `microtrait_trait-name3` = readr::col_character()
                                                          )
                                )
# allsubstrates = substrate2rule %>% pull(`microtrait_substrate-name`)
# 5. read rule to trait
rule2trait_temp = readr::read_delim("./data-raw/microtrait_rule2trait.txt",
                               delim = "\t",
                               col_names = T,
                               col_types = readr::cols(`microtrait_rule-name` = readr::col_factor(ordered = FALSE),
                                                       `microtrait_rule-type` = readr::col_character(),
                                                       `microtrait_rule-substrate` = readr::col_character(),
                                                       `microtrait_trait-name1` = readr::col_character(),
                                                       `microtrait_trait-name2` = readr::col_character(),
                                                       `microtrait_trait-name3` = readr::col_character(),
                                                       `microtrait_trait-version` = readr::col_character()
                                                      )
                                  )
rule2trait_temp = rule2trait_temp %>%
  dplyr::inner_join(rules, by = c("microtrait_rule-name" = "microtrait_rule-name")) %>%
  dplyr::select(c("microtrait_rule-name", "microtrait_rule-boolean", "microtrait_rule-type", "microtrait_rule-substrate",
                  "microtrait_trait-name1", "microtrait_trait-name2", "microtrait_trait-name3",
                  "microtrait_trait-version")) %>%
  dplyr::filter(`microtrait_trait-version` == "production")
## expand count_by_substrate rules
rule2trait_countbysubstrate = rule2trait_temp %>%
  dplyr::filter(`microtrait_rule-type` == "count_by_substrate") %>%
  tidyr::separate_rows(`microtrait_rule-substrate`, sep = ";") %>%
  dplyr::left_join(substrate2rule, by = c("microtrait_rule-substrate" = "microtrait_substrate-name"), suffix = c(".rule2trait", ".substrate2rule")) %>%
  dplyr::select(c("microtrait_rule-name", "microtrait_rule-boolean", "microtrait_rule-type", "microtrait_rule-substrate",
                  "microtrait_trait-name1.substrate2rule", "microtrait_trait-name2.substrate2rule", "microtrait_trait-name3.substrate2rule",
                  "microtrait_trait-version")) %>%
  dplyr::rename(`microtrait_trait-name1` = `microtrait_trait-name1.substrate2rule`,
                `microtrait_trait-name2` = `microtrait_trait-name2.substrate2rule`,
                `microtrait_trait-name3` = `microtrait_trait-name3.substrate2rule`)
# rule2trait_countbysubstrate_allsubstrates = rule2trait_countbysubstrate %>% pull(`microtrait_rule-substrate`) %>% unique

## count and binary rules
rule2trait_remain = rule2trait_temp %>%
  dplyr::filter(`microtrait_rule-type` %in% c("count", "binary"))

rule2trait = dplyr::bind_rows(rule2trait_countbysubstrate, rule2trait_remain)

# check, no NAs in trait-names
# rule2trait%>%pull(`microtrait_trait-name1`)%>%is.na()%>%table()
# rule2trait%>%pull(`microtrait_trait-name2`)%>%is.na()%>%table()
# rule2trait%>%pull(`microtrait_trait-name3`)%>%is.na()%>%table()

rule2trait = rule2trait %>%
  dplyr::mutate(`microtrait_trait-name1` = readr::parse_factor(`microtrait_trait-name1`,
                                                               ordered = T,
                                                               levels = traits_listbygranularity[[1]] %>%
                                                                 dplyr::arrange(`microtrait_trait-displayorder`) %>%
                                                                 dplyr::pull(`microtrait_trait-name`) %>% as.character()
                                                              ),
                `microtrait_trait-name2` = readr::parse_factor(`microtrait_trait-name2`,
                                                               ordered = T,
                                                               levels = traits_listbygranularity[[2]] %>%
                                                                 dplyr::arrange(`microtrait_trait-displayorder`) %>%
                                                                 dplyr::pull(`microtrait_trait-name`) %>% as.character()
                                                              ),
                `microtrait_trait-name3` = readr::parse_factor(`microtrait_trait-name3`,
                                                               ordered = T,
                                                               levels = traits_listbygranularity[[3]] %>%
                                                                 dplyr::arrange(`microtrait_trait-displayorder`) %>%
                                                                 dplyr::pull(`microtrait_trait-name`) %>% as.character()
                                                              )
                )

# 1,701
# 6. read performance metrics for the microtraithmms
hmm_performance = readr::read_delim("./data-raw/microtrait_hmm-performance.txt",
                                              delim = "\t",
                                              col_names = T,
                                              col_types = readr::cols(`microtrait_hmm-dbxref_kegg` = readr::col_character(),
                                                                      `microtrait_hmm-dbxref_kegg_score` = readr::col_double(),
                                                                      `microtrait_hmm-dbxref_kegg_fscore` = readr::col_double(),
                                                                      `microtrait_hmm-dbxref_kegg_tpr` = readr::col_double(),
                                                                      `microtrait_hmm-dbxref_kegg_tnr` = readr::col_double(),
                                                                      `microtrait_hmm-dbxref_kegg_fpr` = readr::col_double(),
                                                                      `microtrait_hmm-dbxref_kegg_fnr` = readr::col_double(),
                                                                      `microtrait_hmm-dbxref_kegg_accuracy` = readr::col_double()
                                                                      )
                                             )
hmms_fromrules_wkegg = hmm_names_fromrules %>%
  dplyr::inner_join(hmm_performance, by = c("microtrait_hmm-dbxref_kegg" = "microtrait_hmm-dbxref_kegg"))
hmms_fromrules_wokegg = hmm_names_fromrules %>%
  dplyr::inner_join(hmm_performance, by = c("microtrait_hmm-name" = "microtrait_hmm-dbxref_kegg"))
hmms_fromrules_cazy = hmm_names_fromrules %>%
  dplyr::inner_join(hmms_cazy, by = c("microtrait_hmm-name" = "microtrait_hmm-name")) %>%
  dplyr::select(c("microtrait_hmm-name", "microtrait_hmm-description.x")) %>%
  dplyr::rename(`microtrait_hmm-description` = `microtrait_hmm-description.x`) %>%
  tibble::add_column(`microtrait_hmm-dbxref_kegg`=NA,
                     `microtrait_hmm-dbxref_ec`=NA,
                     `microtrait_hmm-dbxref_TC`=NA,
                     `microtrait_hmm-notes`=NA,
                     `microtrait_hmm-dbxref_kegg_score`=NA,
                     `microtrait_hmm-dbxref_kegg_fscore`=NA,
                     `microtrait_hmm-dbxref_kegg_tpr`=NA,
                     `microtrait_hmm-dbxref_kegg_tnr`=NA,
                     `microtrait_hmm-dbxref_kegg_fpr`=NA,
                     `microtrait_hmm-dbxref_kegg_fnr`=NA,
                     `microtrait_hmm-dbxref_kegg_accuracy`=NA)

hmms_fromrules = bind_rows(hmms_fromrules_wkegg,
                           hmms_fromrules_wokegg,
                           hmms_fromrules_cazy) %>%
  dplyr::mutate(`microtrait_hmm-name` = readr::parse_factor(`microtrait_hmm-name`,
                                                            ordered = T))
#
# a=c(hmms_fromrules_wkegg %>% pull(`microtrait_hmm-dbxref_kegg`) %>% unique(),
#     hmms_fromrules_wokegg %>% pull(`microtrait_hmm-name`) %>% unique())
# b=hmm_performance %>% pull(`microtrait_hmm-dbxref_kegg`)
# should be empty, otherwise there are genes without performance metrics
# setdiff(a,b)
#hmms_fromrules %>%
#  anti_join(hmm_performance, by = c("microtrait_hmm-dbxref_kegg" = "microtrait_hmm-dbxref_kegg"))

#usethis::use_data(hmms, overwrite = TRUE)
#usethis::use_data(substrate2rule, overwrite = TRUE)
usethis::use_data(hmms_fromrules, overwrite = TRUE)
usethis::use_data(rules, overwrite = TRUE)
usethis::use_data(rule2trait, overwrite = TRUE)
usethis::use_data(traits_listbygranularity, overwrite = TRUE)

write.table(rules,
            file = "./data-raw/microtrait_ruleunwrapped.txt",
            sep = "\t", quote = F, row.names = F, col.names = F)
write.table(hmms_fromrules %>% as.data.frame,
            file = "./data-raw/microtrait_hmmsfromrules.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(hmm_names_fromrules %>% as.data.frame,
            file = "./data-raw/microtrait_hmm_names_fromrules.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)

# build shell commands for "rename_addtc2hmm.sh"
rename_addtc2hmm_path = "$bindir/rename_addtc2hmm.sh"
hmm_indir = "$rawhmmdir"
hmm_outdir = "$renamedwTCdir"
#
rename_addtc2hmm_commands = hmms_fromrules %>%
  dplyr::mutate(`microtrait_hmm-descriptionwithdbxref` = paste0(`microtrait_hmm-description`,
                                                     "//dbxref_kegg", `microtrait_hmm-dbxref_kegg`,
                                                     "//dbxref_ec", `microtrait_hmm-dbxref_ec`,
                                                     "//dbxref_TC", `microtrait_hmm-dbxref_TC`
  )) %>%
  dplyr::mutate(rename_addtc2hmm_command = paste0(rename_addtc2hmm_path, " ",
                                                  hmm_indir, " ",
                                                  hmm_outdir, " ",
                                                  `microtrait_hmm-dbxref_kegg`, " ",
                                                  `microtrait_hmm-name`, " ",
                                                  `microtrait_hmm-dbxref_kegg_score`, " ",
                                                  "\"", `microtrait_hmm-descriptionwithdbxref`, "\"")
  ) %>%
  pull(`rename_addtc2hmm_command`)
write.table(rename_addtc2hmm_commands,
            file = "./data-raw/rename_addtc2hmm_commands.sh", sep = "\t", quote = F, row.names = F, col.names = F)
