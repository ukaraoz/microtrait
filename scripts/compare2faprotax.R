library(dplyr)
library(readr)
library(janitor)

base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))
hmm_matrix = genomeset_results_wmetadata[["hmm_matrix"]]
rule_matrix = genomeset_results_wmetadata[["rule_matrix"]]
hmmandrule_matrix = hmm_matrix %>% dplyr::select(`id`:`OGT`) %>% dplyr::left_join(rule_matrix, by = c("id" = "id"))

source("/Users/ukaraoz/Work/microtrait/code/microtrait/scripts/trait_check.R")
faprotax_traits = readr::read_tsv(file = "./data-raw/faprotax_traits.txt", comment='#') %>%
  dplyr::filter(select == 1) %>%
  dplyr::pull(`faprotax_trait`)

faprotax_traits_matchstats = list()
for(i in 1:length(faprotax_traits)) {
  temp = faprotax_compare_helper(faprotax_traits[i])
  faprotax_traits_matchstats[[faprotax_traits[i]]] = data.frame(faprotax_trait = faprotax_traits[i],
                                               faprotax_taxa = temp[["ntaxa"]],
                                               faprotax_undetermined = temp[["ranks"]]["undetermined"],
                                               faprotax_species = temp[["ranks"]]["species"],
                                               faprotax_genus = temp[["ranks"]]["genus"],
                                               faprotax_family = temp[["ranks"]]["family"],
                                               faprotax_order = temp[["ranks"]]["order"],
                                               faprotax_class = temp[["ranks"]]["class"],
                                               faprotax_phylum = temp[["ranks"]]["phylum"],
                                               NCBI_uniquetaxa = temp[["ncbimatched_ntaxa"]],
                                               NCBI_species = temp[["ncbimatched_ranks"]]["species"],
                                               NCBI_genus = temp[["ncbimatched_ranks"]]["genus"],
                                               NCBI_family = temp[["ncbimatched_ranks"]]["family"],
                                               NCBI_order = temp[["ncbimatched_ranks"]]["order"],
                                               NCBI_class = temp[["ncbimatched_ranks"]]["class"],
                                               NCBI_phylum = temp[["ncbimatched_ranks"]]["phylum"],
                                               IMGgenomes_uniquetaxa = temp[["ncbimatched_genomes_ntaxa"]],
                                               IMGgenomes_species = temp[["ncbimatched_genomes_ranks"]]["species"],
                                               IMGgenomes_genus = temp[["ncbimatched_genomes_ranks"]]["genus"],
                                               IMGgenomes_family = temp[["ncbimatched_genomes_ranks"]]["family"],
                                               IMGgenomes_order = temp[["ncbimatched_genomes_ranks"]]["order"],
                                               IMGgenomes_class = temp[["ncbimatched_genomes_ranks"]]["class"],
                                               IMGgenomes_phylum = temp[["ncbimatched_genomes_ranks"]]["phylum"])
}
faprotax_traits_matchstats[["nitrite_denitrification"]][-1] =
  faprotax_traits_matchstats[["nitrite_denitrification"]][-1] +
  faprotax_traits_matchstats[["nitrate_denitrification"]][-1]
faprotax_traits_matchstats[["nitrous_oxide_denitrification"]][-1] =
  faprotax_traits_matchstats[["nitrous_oxide_denitrification"]][-1] +
  faprotax_traits_matchstats[["nitrite_denitrification"]][-1]
faprotax_traits_matchstats[["denitrification"]][-1] =
  faprotax_traits_matchstats[["denitrification"]][-1] +
  faprotax_traits_matchstats[["nitrate_denitrification"]][-1] +
  faprotax_traits_matchstats[["nitrite_denitrification"]][-1] +
  faprotax_traits_matchstats[["nitrous_oxide_denitrification"]][-1]
faprotax_traits_matchstats[["nitrite_ammonification"]][-1] =
  faprotax_traits_matchstats[["nitrite_ammonification"]][-1] +
  faprotax_traits_matchstats[["nitrate_ammonification"]][-1]
faprotax_traits_matchstats[["nitrate_respiration"]][-1] =
  faprotax_traits_matchstats[["nitrate_respiration"]][-1] +
  faprotax_traits_matchstats[["nitrate_ammonification"]][-1] +
  faprotax_traits_matchstats[["nitrate_denitrification"]][-1]
faprotax_traits_matchstats[["nitrite_respiration"]][-1] =
  faprotax_traits_matchstats[["nitrite_respiration"]][-1] +
  faprotax_traits_matchstats[["nitrate_ammonification"]][-1] +
  faprotax_traits_matchstats[["nitrite_ammonification"]][-1] +
  faprotax_traits_matchstats[["anammox"]][-1]
faprotax_traits_matchstats[["dark_oxidation_of_sulfur_compounds"]][-1] =
  faprotax_traits_matchstats[["dark_oxidation_of_sulfur_compounds"]][-1] +
  faprotax_traits_matchstats[["dark_sulfide_oxidation"]][-1] +
  faprotax_traits_matchstats[["dark_sulfur_oxidation"]][-1] +
  faprotax_traits_matchstats[["dark_thiosulfate_oxidation"]][-1]
faprotax_traits_matchstats[["hydrogenotrophic_methanogenesis"]][-1] =
  faprotax_traits_matchstats[["hydrogenotrophic_methanogenesis"]][-1] +
  faprotax_traits_matchstats[["methanogenesis_by_CO2_reduction_with_H2"]][-1]
#  #faprotax_traits_matchstats[["methanogenesis_by_reduction_of_methyl_compounds_with_H2"]][-1]
faprotax_traits_matchstats[["methanogenesis"]][-1] =
  faprotax_traits_matchstats[["methanogenesis"]][-1] +
  faprotax_traits_matchstats[["hydrogenotrophic_methanogenesis"]][-1] +
  faprotax_traits_matchstats[["methanogenesis_by_disproportionation_of_methyl_groups"]][-1] +
  faprotax_traits_matchstats[["acetoclastic_methanogenesis"]][-1]
  #faprotax_traits_matchstats[["methanogenesis_using_formate"]][-1] +
  #faprotax_traits_matchstats[["methanogenesis_by_CO2_reduction_with_H2"]][-1]
  #faprotax_traits_matchstats[["methanogenesis_by_reduction_of_methyl_compounds_with_H2"]][-1]
faprotax_traits_matchstats[["methylotrophy"]][-1] =
  faprotax_traits_matchstats[["methylotrophy"]][-1] +
  faprotax_traits_matchstats[["methanol_oxidation"]][-1] +
  faprotax_traits_matchstats[["methanotrophy"]][-1]
faprotax_traits_matchstats[["anoxygenic_photoautotrophy"]][-1] =
  faprotax_traits_matchstats[["anoxygenic_photoautotrophy"]][-1] +
  #faprotax_traits_matchstats[["anoxygenic_photoautotrophy_H2_oxidizing"]][-1] +
  faprotax_traits_matchstats[["anoxygenic_photoautotrophy_S_oxidizing"]][-1] +
  faprotax_traits_matchstats[["anoxygenic_photoautotrophy_Fe_oxidizing"]][-1]
faprotax_traits_matchstats[["photoheterotrophy"]][-1] =
  faprotax_traits_matchstats[["photoheterotrophy"]][-1] +
  faprotax_traits_matchstats[["aerobic_anoxygenic_phototrophy"]][-1]

faprotax_traits_matchstats = do.call(rbind, faprotax_traits_matchstats)
write.table(faprotax_traits_matchstats,
            file = "/Users/ukaraoz/Work/microtrait/code/microtrait/data-raw/faprotax_traits_matchstats.xls",
            row.names = F, col.names = T, quote = F, sep = "\t")

trait_check_list = trait_check_list()
genes = c(as.character(unlist(trait_check_list[["gene_check"]][["nitrogen"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["sulfur"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["methanogenesis"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["methanotrophy"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["phototrophy"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["sulfur"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["iron"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["fumarate_respiration"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["arsenic"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["selenium"]])),
          as.character(unlist(trait_check_list[["gene_check"]][["carbonfixation"]])))
rules = as.character(unlist(trait_check_list$rule_check))


faprotax_traits = setdiff(faprotax_traits, "photoheterotrophy")
faprotax_wtrait_ncbimatched_genomes = assert_faprotaxtraits_forgenomes(faprotax_traits)
faprotax_wtrait_ncbimatched_genomes_sums = faprotax_wtrait_ncbimatched_genomes %>% summarise(across(faprotax_traits, ~ sum(.)))
faprotax_traitcheck = hmmandrule_matrix %>%
  dplyr::left_join(faprotax_wtrait_ncbimatched_genomes, by = c("id" = "id")) %>%
  dplyr::mutate_at(vars(faprotax_traits), funs(replace(., which(is.na(.)), 0))) %>%
  dplyr::rename(mtrD = `mtrD.x`) %>%
  dplyr::select(c("id", "NCBI_Organism Name", "NCBI_Superkingdom", "NCBI_Phylum", "NCBI_Class",
                  "NCBI_Order", "NCBI_Family", "NCBI_Genus", "NCBI_Species",
                  faprotax_traits, rules, genes))
write.table(faprotax_traitcheck,
            file = "/Users/ukaraoz/Work/microtrait/code/microtrait/data-raw/faprotax_traitcheck.xls",
            quote = F, row.names = F, col.names = T, sep = "\t")

#faprotax_traits1 = c("nitrate_denitrification", "nitrite_denitrification", "nitrous_oxide_denitrification", "nitrate_ammonification",
#                    "nitrite_ammonification", "nitrogen_fixation", "sulfate_respiration", "sulfite_respiration",
#                    "sulfur_respiration", "thiosulfate_respiration", "dark_sulfide_oxidation", "dark_sulfite_oxidation",
#                    "dark_sulfur_oxidation", "dark_thiosulfate_oxidation", "hydrogenotrophic_methanogenesis", "methanogenesis_by_CO2_reduction_with_H2",
#                    "methanogenesis_by_disproportionation_of_methyl_groups", "acetoclastic_methanogenesis", "oxygenic_photoautotrophy",
#                    "anoxygenic_photoautotrophy_S_oxidizing", "anoxygenic_photoautotrophy_Fe_oxidizing", "aerobic_anoxygenic_phototrophy",
#                    "iron_respiration", "dark_iron_oxidation", "chlorate_reducers", "fumarate_respiration", "reductive_acetogenesis")

microtrait_functionalgroups_forfaprotax = define_functionalgroups_forfaprotax(hmmandrule_matrix, faprotax_traits)
hits_table = microtrait_functionalgroups_forfaprotax %>%
  dplyr::left_join(faprotax_wtrait_ncbimatched_genomes, by = c("id" = "id")) %>%
  dplyr::mutate_each(funs(replace(., which(is.na(.)), 0)))

write.table(hits_table,
            file = "/Users/ukaraoz/Work/microtrait/code/microtrait/data-raw/faprotax_microtrait_comparison_hits_table.xls",
            row.names = F, col.names = T, quote = F, sep = "\t")

results = lapply(1:length(faprotax_traits),
                 function(i) {returnList = crosstable(hits_table, faprotax_traits[i])[[1]]
                              returnList
                 })
results = do.call(rbind, results)
write.table(results,
            file = "/Users/ukaraoz/Work/microtrait/code/microtrait/data-raw/faprotax_microtrait_comparison.xls",
           row.names = F, col.names = T, quote = F, sep = "\t")



assert_faprotaxtraits_forgenomes <- function(faprotax_traits) {
  faprotax_wtrait_ncbimatched_genomes = faprotax %>%
    dplyr::filter(`functional group` %in% faprotax_traits) %>%
    dplyr::filter(`ncbi_taxa` != "unmatched") %>%
    dplyr::inner_join(faprotax2genomes, by = c("taxa" = "faprotax_taxa")) %>%
    #dplyr::filter(`id` == "2505119042" & `functional group` == "sulfate_respiration")
    dplyr::select(c("id", `functional group`)) %>%
    dplyr::distinct(`id`, `functional group`) %>%
    dplyr::mutate(present = as.integer(1)) %>%
    tidyr::pivot_wider(names_from = "functional group",
                       values_from = "present",
                       values_fill = 0) %>%
    dplyr::select(c("id", faprotax_traits))

  faprotax_wtrait_ncbimatched_genomes = faprotax_wtrait_ncbimatched_genomes %>%
    dplyr::mutate(`nitrite_denitrification` =
                    case_when(`nitrite_denitrification`==1~as.integer(1),
                              `nitrate_denitrification`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`nitrous_oxide_denitrification` =
                    case_when(`nitrous_oxide_denitrification`==1~as.integer(1),
                              `nitrite_denitrification`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`denitrification` =
                    case_when(`denitrification`==1~as.integer(1),
                              `nitrate_denitrification`==1~as.integer(1),
                              `nitrite_denitrification`==1~as.integer(1),
                              `nitrous_oxide_denitrification`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`nitrite_ammonification` =
                    case_when(`nitrite_ammonification`==1~as.integer(1),
                              `nitrate_ammonification`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`nitrate_respiration` =
                    case_when(`nitrate_respiration`==1~as.integer(1),
                              `nitrate_ammonification`==1~as.integer(1),
                              `nitrate_denitrification`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`nitrite_respiration` =
                    case_when(`nitrite_respiration`==1~as.integer(1),
                              `nitrate_ammonification`==1~as.integer(1),
                              `nitrite_ammonification`==1~as.integer(1),
                              `anammox`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`dark_oxidation_of_sulfur_compounds` =
                    case_when(`dark_oxidation_of_sulfur_compounds`==1~as.integer(1),
                              `dark_sulfide_oxidation`==1~as.integer(1),
                              `dark_sulfur_oxidation`==1~as.integer(1),
                              `dark_thiosulfate_oxidation`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`hydrogenotrophic_methanogenesis` =
                    case_when(`hydrogenotrophic_methanogenesis`==1~as.integer(1),
                              `methanogenesis_by_CO2_reduction_with_H2`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`methanogenesis` =
                    case_when(`methanogenesis`==1~as.integer(1),
                              `hydrogenotrophic_methanogenesis`==1~as.integer(1),
                              `methanogenesis_by_disproportionation_of_methyl_groups`==1~as.integer(1),
                              `acetoclastic_methanogenesis`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`methylotrophy` =
                    case_when(`methylotrophy`==1~as.integer(1),
                              `methanol_oxidation`==1~as.integer(1),
                              `methanotrophy`==1~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(`anoxygenic_photoautotrophy` =
                    case_when(`anoxygenic_photoautotrophy`==1~as.integer(1),
                              `anoxygenic_photoautotrophy_S_oxidizing`==1~as.integer(1),
                              `anoxygenic_photoautotrophy_Fe_oxidizing`==1~as.integer(1),
                              TRUE~as.integer(0))) # %>%
    #dplyr::mutate(`photoheterotrophy` =
    #                case_when(`photoheterotrophy`==1~as.integer(1),
    #                          `aerobic_anoxygenic_phototrophy`==1~as.integer(1),
    #                          TRUE~as.integer(0)))
    return(faprotax_wtrait_ncbimatched_genomes)
}

crosstable <- function(hits_table, faprotax_trait) {
  temp = hits_table %>%
    dplyr::select(c("id", faprotax_trait, paste0("microtrait_", faprotax_trait))) %>%
    dplyr::mutate(faprotax = get(faprotax_trait),
                  microtrait = get(paste0("microtrait_", faprotax_trait))) %>%
    dplyr::select(c("id", "faprotax", "microtrait"))
  crosstab = temp %>% janitor::tabyl(faprotax, microtrait)

  #adorn_totals(c("row", "col")) %>%
  #  adorn_percentages("row") %>%
  #  adorn_pct_formatting(rounding = "half up", digits = 0) %>%
  #  adorn_ns() %>%
  #  adorn_title("combined")
  performance = list()
  performance$faprotax_trait = faprotax_trait
  performance$TN = crosstab[[2]][1]
  performance$FN = crosstab[[2]][2]
  if(length(crosstab) == 3) {
    performance$FP = crosstab[[3]][1]
    performance$TP = crosstab[[3]][2]
  } else {
    performance$FP = 0
    performance$TP = 0
  }
  performance$P = performance$TP + performance$FN
  performance$N = performance$TN + performance$FP
  performance$TPR = performance$TP/performance$P*100
  performance$TNR = performance$TN/performance$N*100

  ids = list()
  ids$TP = temp %>% dplyr::filter(faprotax==1&microtrait==1) %>% dplyr::pull(id)
  ids$FP = temp %>% dplyr::filter(faprotax==1&microtrait==0) %>% dplyr::pull(id)
  ids$TN = temp %>% dplyr::filter(faprotax==0&microtrait==0) %>% dplyr::pull(id)
  ids$FN = temp %>% dplyr::filter(faprotax==0&microtrait==1) %>% dplyr::pull(id)
  result = list(performance = performance, ids = ids)
  return(result)
}

define_functionalgroups_forfaprotax <- function(hmmandrule_matrix, faprotax_traits) {
  faprotax_groups = hmmandrule_matrix %>%
    dplyr::mutate(microtrait_denitrification =
                    #case_when(((narG==1|narH==1|narI_dsrM==1|napA==1|napB==1))~as.integer(1),
                    case_when(((`narGHI_and`==1|`napAB_and`==1|`nirSK_or`==1|`nosZ`==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_nitrate_denitrification =
                    #case_when(((narG==1|narH==1|narI_dsrM==1|napA==1|napB==1))~as.integer(1),
                    case_when(((`narGHI_and`==1|`napAB_and`==1))~as.integer(1),
                               TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_nitrite_denitrification =
                    #case_when(((narG==1|narH==1|narI_dsrM==1|napA==1|napB==1|nirS==1|nirK==1))~as.integer(1),
                    case_when(((`narGHI_and`==1|`napAB_and`==1|`nirSK_or`==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_nitrous_oxide_denitrification =
                    #case_when(((narG==1|narH==1|narI_dsrM==1|napA==1|napB==1|nirS==1|nirK==1|nosZ==1))~as.integer(1),
                    case_when(((`narGHI_and`==1|`napAB_and`==1|`nirSK_or`==1|`nosZ`==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_nitrate_respiration =
                    #case_when(((narG==1|narH==1|narI_dsrM==1|napA==1|napB==1))~as.integer(1),
                    case_when(((`narGHI_or`==1|`napAB_or`==1)|(nirB==1|nirD==1|nrfA==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_nitrate_ammonification =
                    case_when((((`narGHI_or`==1|`napAB_or`==1)&(nirB==1|nirD==1|nrfA==1)))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_nitrite_ammonification =
                    case_when((nirB==1|nirD==1|nrfA==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_nitrite_respiration =
                    case_when((nirB==1|nirD==1|nrfA==1|`nirSK_or`==1|hzsA==1|hzsB==1|hzsC==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_nitrogen_fixation =
                    case_when(((nifH==1|vnfH==1|nifK==1|nifD==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_aerobic_ammonia_oxidation =
                    case_when((`pmoA-amoA`==1|`pmoB-amoB`==1|`pmoC-amoC`==1|hao==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_aerobic_nitrite_oxidation =
                    case_when((narG==1|narH==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_anammox =
                    case_when((hzsA==1|hzsB==1|hzsC==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_sulfate_respiration =
                    case_when(sat==1~as.integer(1), TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_sulfite_respiration =
                    case_when(((dsrA==1|dsrB==1|dsrE==1|dsrF==1|dsrH==1|asrA==1|asrB==1|asrC==1|mccA==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_sulfur_respiration =
                    case_when((sat==1|aprA==1|aprB==1|sreA==1|sreB==1|sreC==1|sqr==1|hydG==1|hydB_4==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_thiosulfate_respiration =
                    case_when((`phsABC_or`==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_dark_oxidation_of_sulfur_compounds =
                    case_when(((fccA==1|fccB==1)|(dsrA==1|dsrB==1|sor==1)|(soxB==1|soxC==1|soxD==1|tsdA==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_dark_sulfide_oxidation =
                    case_when(((fccA==1|fccB==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_dark_sulfite_oxidation =
                    case_when((sorA==1|sorB==1|soeA==1|soeB==1|soeC==1|aprA==1|aprB==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_dark_sulfur_oxidation =
                    case_when(((dsrA==1|dsrB==1|sor==1|soxB==1|soxC==1|soxD==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_dark_thiosulfate_oxidation =
                    case_when((soxB==1|soxC==1|soxD==1|tsdA==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_methanogenesis =
                  case_when(((`methylCoM+CoB->CoM+CH4_or`==1))~as.integer(1),
                            TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_hydrogenotrophic_methanogenesis =
                  #case_when(((mcrA==1&mcrB==1&mcrG==1&mcr2==1&mcrC==1&mcrD==1&fwdA==1&fwdB==1&fwdC==1&fwdD==1&fwdF==1&fwdH==1&fwdG==1&fwdE==1&ftr==1&mtdch==1&hmd==1&mer==1))~as.integer(1),
                  case_when(((`MF+CO2->formylMF_or`==1&
                              `formylMF+H4MPT->MF+formylH4MPT`==1&
                              `formylH4MPT->methenylH4MPT`==1&
                              `methyleneH4MPT->methyl-H4SPT`==1&
                              `methylCoM+CoB->CoM+CH4_or`==1))~as.integer(1),
                                 TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_methanogenesis_by_CO2_reduction_with_H2 =
                    case_when(((`methylCoM+CoB->CoM+CH4_and`==1&ftr==1&mtdch==1&mer==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_methanogenesis_by_disproportionation_of_methyl_groups =
                    case_when(((`methylCoM+CoB->CoM+CH4_and`==1)&(mtbA==1|mtmB==1|mtmC==1|mtbB==1|mtbC==1|mttB==1|mttC==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_acetoclastic_methanogenesis =
                    case_when(((mcrA==1&mcrB==1&mcrG==1&mcr2==1&mcrC==1&mcrD==1&
                                  ackA==1&cdhC==1&acsC==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_methylotrophy =
                    case_when((`CH4->CH3OH`==1|`CH3OH->CH2O`==1)~as.integer(1),
                            TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_methanol_oxidation =
                    case_when((`CH3OH->CH2O`==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_methanotrophy =
                    case_when((`CH4->CH3OH`==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_reductive_acetogenesis =
                    case_when(((ftfL==1&folD==1&metF==1&acsE==1&acsB==1&acsC==1&acsD==1&ackA==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_oxygenic_photoautotrophy =
                    case_when((`oxygenic photosystem I core`==1&`oxygenic photosystem II core`==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_anoxygenic_photoautotrophy =
                    case_when(((`anoxygenic photosystem I RC`==1|
                                  `anoxygenic photosystem II RC`==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_anoxygenic_photoautotrophy_S_oxidizing =
                    case_when(((`anoxygenic photosystem I RC`==1|
                                  `anoxygenic photosystem II RC`==1)&
                                 (`sulfite->sulfate-oxidation`==1|
                                    `sulfide->sulfur-oxidation`==1|
                                    `sulfur->sulfite-oxidation`==1|
                                    `sulfur->sulfate-oxidation`==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_anoxygenic_photoautotrophy_Fe_oxidizing =
                    case_when(((`anoxygenic photosystem I RC`==1|
                                  `anoxygenic photosystem II RC`==1)&
                                 (pioA==1|pioB==1|pioC==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%         #photoheterotrophy:  aerobic_anoxygenic_phototrophy++
    dplyr::mutate(microtrait_aerobic_anoxygenic_phototrophy =
                    case_when(((bchG==1&bciB==0&por==0&bciC==0&bchV==0&bchU==0&bchK==0&bciD==0))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_iron_respiration =
                    case_when((mtrA==1|mtrB==1|mtrC==1|mtrD.x==1|mtrE==1|mtrF==1|omcA==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_dark_iron_oxidation =
                    case_when((pioA==1|pioB==1|pioC==1|foxE==1|foxY==1|foxZ==1|cyc2==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_chlorate_reducers =
                    case_when((serA==1|serB==1|serC==1)~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_fumarate_respiration =
                    case_when(((`fumarate->succinate-frd_or`==1|
                                `fumarate->succinate-frdNADH_or`==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_fermentation =
                    case_when(((`acetoacetyl-CoA->acetone`==1|
                                `pyruvate->lactate and ethanol`==1|
                                `pyruvate->alcohols`==1|
                                `pyruvate->SCFAs`==1|
                                `pyruvate->butyleneglycol`==1|
                                `pyruvate->butanol`==1|
                                `acetylene->acetate`==1|
                                `acetylene->acetate+ethanol`==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::mutate(microtrait_aerobic_chemoheterotrophy =
                    case_when(((`nuo_and`==1|
                                `succinate->fumarate_key`==1|
                                `cytochromecreductase_or`==1|
                                `cytochromeoubiquinol_and`==1|
                                `cytochromec_and`==1|
                                `cytochromeaa3600menaquinol_and`==1|
                                `cytochromeccbb3oxidase_and`==1))~as.integer(1),
                              TRUE~as.integer(0))) %>%
    dplyr::select(c(`id`,
                    "NCBI_Organism Name", "NCBI_Superkingdom", "NCBI_Phylum", "NCBI_Class",
                    "NCBI_Order", "NCBI_Family", "NCBI_Genus", "NCBI_Species",
                    "Ecosystem", "Ecosystem_Category", "Ecosystem_Type",
                    "Ecosystem_Subtype","Specific Ecosystem",
                    c(rbind(faprotax_traits, paste0("microtrait_", faprotax_traits, sep = "")))
                    #paste0("microtrait_", faprotax_traits))
                  ))
  return(faprotax_groups)
}

faprotax_compare_helper <- function(faprotax_trait, genes = NULL) {
  faprotax_wtrait = faprotax %>% dplyr::filter(`functional group` == faprotax_trait)

  faprotax_wtrait_taxa_ranks = faprotax_wtrait %>% dplyr::select("taxalevel") %>% table()
  faprotax_wtrait_taxa_nrank = rep(NA, 7)
  names(faprotax_wtrait_taxa_nrank) = c("undetermined", "species", "genus", "family", "order", "class", "phylum")
  faprotax_wtrait_taxa_nrank[names(faprotax_wtrait_taxa_ranks)] = faprotax_wtrait_taxa_ranks
  faprotax_wtrait_taxa_nrank[is.na(faprotax_wtrait_taxa_nrank)] = 0

  faprotax_wtrait_ncbimatched = faprotax_wtrait %>% dplyr::filter(`ncbi_taxa` != "unmatched")
  faprotax_wtrait_ncbimatched_ranks = faprotax_wtrait_ncbimatched %>% dplyr::select(c("rank")) %>% table()
  faprotax_wtrait_ncbimatched_nrank = rep(NA, 6)
  names(faprotax_wtrait_ncbimatched_nrank) = c("species", "genus", "family", "order", "class", "phylum")
  faprotax_wtrait_ncbimatched_nrank[names(faprotax_wtrait_ncbimatched_ranks)] = faprotax_wtrait_ncbimatched_ranks
  faprotax_wtrait_ncbimatched_nrank[is.na(faprotax_wtrait_ncbimatched_nrank)] = 0

  faprotax_wtrait_ncbimatched_genomes = faprotax_wtrait %>%
    dplyr::filter(`ncbi_taxa` != "unmatched") %>%
    dplyr::inner_join(faprotax2genomes, by = c("taxa" = "faprotax_taxa"))
  faprotax_wtrait_ncbimatched_genomes_ranks = faprotax_wtrait_ncbimatched_genomes %>% dplyr::select("rank.x") %>% table()
  faprotax_wtrait_ncbimatched_genomes_nrank = rep(NA, 6)
  names(faprotax_wtrait_ncbimatched_genomes_nrank) = c("species", "genus", "family", "order", "class", "phylum")
  faprotax_wtrait_ncbimatched_genomes_nrank[names(faprotax_wtrait_ncbimatched_genomes_ranks)] = faprotax_wtrait_ncbimatched_genomes_ranks
  faprotax_wtrait_ncbimatched_genomes_nrank[is.na(faprotax_wtrait_ncbimatched_genomes_nrank)] = 0

  faprotax_wtrait_ncbimatched_genomes_wgenes = faprotax_wtrait_ncbimatched_genomes %>%
    dplyr::inner_join(hmm_matrix, by = c("id" = "id")) %>%
    dplyr::select(c("taxa", "taxalevel", "ncbi_taxa.x", "Organism Name", "id", genes))

  faprotax_wtrait_ncbimatched_genomes_wgenes_collapsed = faprotax_wtrait_ncbimatched_genomes_wgenes %>%
    group_by(`taxa`) %>%
    mutate_at(genes, funs(as.numeric(as.character(.)))) %>%
    summarise(across(genes, ~ sum(.))) %>%
    mutate_at(genes, funs(case_when(.>0 ~ 1,
                                    .<=0 ~ 0)))
  result = list()
  result[["faprotax_trait"]] = faprotax_trait
  result[["ntaxa"]] = sum(faprotax_wtrait_taxa_ranks)
  result[["ranks"]] = faprotax_wtrait_taxa_nrank
  # number of faprotax taxa after ncbimatch
  result[["ncbimatched_ntaxa"]] = faprotax_wtrait_ncbimatched %>% dplyr::select(`taxa`) %>% n_distinct()
  result[["ncbimatched_ranks"]] = faprotax_wtrait_ncbimatched_nrank
  # number of faprotax taxa after ncbimatch and genome match
  result[["ncbimatched_genomes_ntaxa"]] = faprotax_wtrait_ncbimatched_genomes %>% dplyr::select(`taxa`) %>% n_distinct()
  result[["ncbimatched_genomes_ranks"]] = faprotax_wtrait_ncbimatched_genomes_nrank
  result[["ncbimatched_genomes_wgenes"]] = faprotax_wtrait_ncbimatched_genomes_wgenes
  result[["ncbimatched_genomes_wgenes_collapsed"]] = faprotax_wtrait_ncbimatched_genomes_wgenes_collapsed

  return(result)
}




#########################

#trait_check_list[["faprotax_check"]][[]]
#results = list()
#for(i in 42:50) {
#  results[[i]] = faprotax_compare(faprotax_traits[i], genes)[["ncbimatched_genomes_wgenes_collapsed"]] %>%
#    dplyr::mutate(faprotax_trait = faprotax_traits[i]) %>%
#    dplyr::select(faprotax_trait, everything())
#  #cat(faprotax_traits[i], "\t", result$ntaxa, "\t", result$ncbimatched_ntaxa, "\t", result$ncbimatched_ngenomes, "\n")
#}
#results = do.call(bind_rows, results)
#write.table(results %>% as.data.frame(),
#            file = "/Users/ukaraoz/Work/microtrait/code/microtrait/data-raw/faprotax_compare_results_arsenic_selenium_iron.xls",
#            row.names = F, col.names = T, quote = F, sep = "\t")

#faprotax_trait = "dark_thiosulfate_oxidation"
#rule = trait_check_list[["faprotax_check"]][[faprotax_trait]]
faprotax_wtrait_ncbimatched_genomes = faprotax_wtrait %>%
  dplyr::filter(`ncbi_taxa` != "unmatched") %>%
  dplyr::inner_join(faprotax2genomes, by = c("taxa" = "faprotax_taxa"))
#genes = as.character(unlist(trait_check_list[["gene_check"]][["nitrogen"]]))
#genes = as.character(unlist(trait_check_list[["gene_check"]][["sulfur"]]))
#genes = c(as.character(unlist(trait_check_list[["gene_check"]][["methanogenesis"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["methanotrophy"]])))
#genes = c(as.character(unlist(trait_check_list[["gene_check"]][["phototrophy"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["sulfur"]][["sulfite->sulfate-oxidation"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["sulfur"]][["sulfide->sulfur-oxidation"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["sulfur"]][["sulfur->sulfite-oxidation"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["sulfur"]][["sulfur->sulfate-oxidation"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["iron"]])))
#genes = c(as.character(unlist(trait_check_list[["gene_check"]][["fumarate_respiration"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["arsenic"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["selenium"]])),
#          as.character(unlist(trait_check_list[["gene_check"]][["iron"]])))
