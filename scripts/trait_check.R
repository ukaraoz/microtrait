#' Return genes that are part of a rule.
#'
#' @return result
#'
#' @export trait_check_list
trait_check_list <- function() {
    rule_check = list()
    rule_check[["nitrogen"]][["denitrification"]][["steps"]][["full"]] = c("nitrate->nitrite_dissimilatory_full",
                                                             "nitrite->nitricoxide_full",
                                                             "nitricoxide->nitrousoxide_full",
                                                             "nitrousoxide->dinitrogen")
    rule_check[["nitrogen"]][["denitrification"]][["steps"]][["partial"]] = c("nitrate->nitrite_dissimilatory_partial",
                                                                "nitrite->nitricoxide_partial",
                                                                "nitricoxide->nitrousoxide_partial",
                                                                "nitrousoxide->dinitrogen")
    rule_check[["nitrogen"]][["denitrification"]][["endproduct"]] = c("nitrate->nitrite_dissimilatory_full",
                                                        "nitrate->nitricoxide_complete",
                                                        "nitrate->nitrousoxide_complete",
                                                        "nitrate->dinitrogen_complete")
    rule_check[["nitrogen"]][["nitrogen fixation"]][["full"]] = "dinitrogen->ammonia_full"
    rule_check[["nitrogen"]][["nitrogen fixation"]][["partial"]] = "dinitrogen->ammonia_partial"
    rule_check[["nitrogen"]][["nitrogen fixation"]][["key"]] = "dinitrogen->ammonia_key"
    rule_check[["nitrogen"]][["ammonia oxidation"]][["full"]] = c("ammonia->hydroxylamine_full", "hydroxylamine->nitrite")
    rule_check[["nitrogen"]][["ammonia oxidation"]][["partial"]] = c("ammonia->hydroxylamine_partial", "hydroxylamine->nitrite")
    rule_check[["nitrogen"]][["nitrite oxidation"]][["full"]] = c("nitrite oxidation_full")
    rule_check[["nitrogen"]][["nitrite oxidation"]][["partial"]] = c("nitrite oxidation_partial")
    rule_check[["methanogenesis"]][["full"]] = c("methylCoM+CoB->CoM+CH4_and")
    rule_check[["methanogenesis"]][["partial"]] = c("methylCoM+CoB->CoM+CH4_or")
    rule_check[["methanogenesis"]][["key"]] = c("methanogenesis_key")

    rule_check[["methanogenesis"]][["acetoclastic"]][["full"]] = c("acetyl-P<->acetate",
                                                                   "acetyl-P->acetyl-CoA",
                                                                   "acetyl-CoA+methylH4MPT->methylH4SPT+CO_full",
                                                                   "methylH4SPT+CoM->methylCoM+methylH4MPT_and",
                                                                   "methylCoM+CoB->CoM+CH4_and")

    rule_check[["methanogenesis"]][["acetoclastic"]][["partial"]] = c("acetyl-P<->acetate",
                                                                      "acetyl-P->acetyl-CoA",
                                                                      "acetyl-CoA+methylH4MPT->methylH4SPT+CO_partial",
                                                                      "methylH4SPT+CoM->methylCoM+methylH4MPT_or",
                                                                      "methylCoM+CoB->CoM+CH4_or")
    rule_check[["methanogenesis"]][["acetoclastic"]][["key"]] = c("acetyl-P<->acetate",
                                                                  "acetyl-P->acetyl-CoA",
                                                                  "acetyl-CoA+methylH4MPT->methylH4SPT+CO_partial",
                                                                  "methylH4SPT+CoM->methylCoM+methylH4MPT_or",
                                                                  "methanogenesis_key")
    rule_check[["methanogenesis"]][["hydrogenotrophic"]][["full"]] = c("MF+CO2->formylMF_and",
                                                                       "formylMF+H4MPT->MF+formylH4MPT",
                                                                       "formylH4MPT->methenylH4MPT",
                                                                       "methenylH4MPT->methyleneH4MPT",
                                                                       "methyleneH4MPT->methyl-H4SPT")
    rule_check[["methanogenesis"]][["hydrogenotrophic"]][["partial"]] = c("MF+CO2->formylMF_or",
                                                                          "formylMF+H4MPT->MF+formylH4MPT",
                                                                          "formylH4MPT->methenylH4MPT",
                                                                          "methenylH4MPT->methyleneH4MPT",
                                                                          "methyleneH4MPT->methyl-H4SPT")
    rule_check[["methanogenesis"]][["methylotrophic"]][["full"]] = c("methanethiol->methylCoM",
                                                                     "methanol->methylCoM_and",
                                                                     "TMAO->trimethylamine",
                                                                     "monomethylamine->methylCoM-direct_and",
                                                                     "dimethylamine->methylCoM-direct_and",
                                                                     "trimethylamine->methylCoM-direct_and")
    rule_check[["methanogenesis"]][["methylotrophic"]][["partial"]] = c("methanethiol->methylCoM",
                                                                     "methanol->methylCoM_or",
                                                                     "TMAO->trimethylamine",
                                                                     "monomethylamine->methylCoM-direct_or",
                                                                     "dimethylamine->methylCoM-direct_or",
                                                                     "trimethylamine->methylCoM-direct_or")
    rule_check[["sulfur"]]
    rule_check[["ETC"]][["ETC_complexI"]][["full"]] = c("nuo_and")
    rule_check[["ETC"]][["ETC_complexII"]][["full"]] = c("succinate->fumarate_full")
    rule_check[["ETC"]][["ETC_complexIII"]][["full"]] = c("cytochromecreductase_or")
    rule_check[["ETC"]][["ETC_complexIV"]][["full"]] = c("cytochromeoubiquinol_and",
                                                         "cytochromec_and",
                                                         "cytochromeaa3600menaquinol_and",
                                                         "cytochromeccbb3oxidase_and")
    rule_check[["fermentation"]] = c("acetoacetyl-CoA->acetone",
                                     "pyruvate->lactate and ethanol",
                                     "pyruvate->alcohols",
                                     "pyruvate->SCFAs",
                                     "pyruvate->butyleneglycol",
                                     "pyruvate->butanol",
                                     "acetylene->acetate",
                                     "acetylene->acetate+ethanol")
    rule_check[["heterotrophy"]] = c("glycolysis-ED_full",
                                     "glycolysis-EMP_full",
                                     "fructose degradation",
                                     "fucose degradation",
                                     "maltose degradation",
                                     "galactose degradation",
                                     "mannose degradation",
                                     "trehalose degradation",
                                     "galacturonate degradation",
                                     "glucuronate degradation",
                                     "glycerol degradation",
                                     "resorcinol degradation")
    gene_check = list()
    gene_check[["nitrogen"]][["denitrification"]] = c("narG", "narH", "narI_dsrM", "napA", "napB", "nirS", "nirK", "norB", "norC", "norV", "norW", "nosZ")
    gene_check[["nitrogen"]][["nitrogen fixation"]] = c("nifD", "nifK", "nifH", "anfG", "vnfD", "vnfG", "vnfK", "vnfH")
    gene_check[["nitrogen"]][["ammonia oxidation"]] = c("pmoA-amoA", "pmoB-amoB", "pmoC-amoC", "hao")
    gene_check[["nitrogen"]][["nitrite oxidation"]] = c("narG", "narH")
    gene_check[["nitrogen"]][["anammox"]] = c("hzsA", "hzsB", "hzsC", "hdh")
    gene_check[["nitrogen"]][["dissimilatory nitrite reduction"]] = c("nirB", "nirD", "nrfA", "nrfH")

    gene_check[["methanogenesis"]][["lowerbranch"]] = c("tmtrA", "tmtrB", "tmtrC", "tmtrD", "tmtrE", "tmtrF", "tmtrG", "tmtrH", "mcrA", "mcrB", "mcrG", "mcr2", "mcrC", "mcrD")
    gene_check[["methanogenesis"]][["acetoclastic"]] = c("ackA", "pta1", "pta2", "cdhC", "acsC", "acsD", gene_check[["methanogenesis"]][["lowerbranch"]])
    gene_check[["methanogenesis"]][["hydrogenotrophic"]] = c("fwdA", "fwdB", "fwdC", "fwdD", "fwdF", "fwdH", "fwdG", "fwdE", "ftr", "mtdch", "hmd", "mer", gene_check[["methanogenesis"]][["lowerbranch"]])
    gene_check[["methanogenesis"]][["methylotrophic"]] = c("mmtsA", "mtaB", "mtaC", "torA", "torZ", "torC", "torY", "mtbA", "mtmB", "mtmC", "mtbA", "mtbB", "mtbC", "mtbA", "mttB", "mttC", gene_check[["methanogenesis"]][["lowerbranch"]])

    gene_check[["methanotrophy"]][["methaneoxidation"]] = c("mmoX", "mmoY", "mmoZ", "mmoC", "mmoD", "pmoA-amoA", "pmoB-amoB", "pmoC-amoC")
    gene_check[["methanotrophy"]][["methanoloxidation"]] = c("mxaF", "mxaJ", "mxaG", "mxaI", "mxaA", "mxaC", "mxaD", "mxaK", "mxaL", "xoxF")

    gene_check[["sulfur"]][["sulfite->sulfide-reduction"]] = c("dsrA","dsrB", "dsrE", "dsrF", "dsrH", "asrA", "asrB", "asrC", "mccA")
    gene_check[["sulfur"]][["sulfate->sulfite-reduction"]] = c("sat", "aprA", "aprB")
    gene_check[["sulfur"]][["sulfur->sulfide-reduction"]] = c("sreA", "sreB", "sreC", "hydG", "hydB_4")
    gene_check[["sulfur"]][["sulfite->sulfate-oxidation"]] = c("aprA", "aprB", "sat", "sorA", "sorB", "soeA", "soeB", "soeC")
    gene_check[["sulfur"]][["sulfide->sulfur-oxidation"]] = c("fccA", "fccB", "sqr")
    gene_check[["sulfur"]][["sulfur->sulfite-oxidation"]] = c("dsrA", "dsrB", "dsrD", "dsrE", "dsrF", "dsrH")
    gene_check[["sulfur"]][["sulfur->sulfate-oxidation"]] = c("dsrA", "dsrB", "dsrD", "dsrE", "dsrF", "dsrH", "aprA", "aprB", "sat")
    gene_check[["sulfur"]][["sulfur->sulfide+sulfite+thiosulfate-disproportionation"]] = c("sor")
    gene_check[["sulfur"]][["thiosulfate->sulfite+sulfide-disproportionation"]] = c("phsA", "phsB", "phsC", "sseA")
    gene_check[["sulfur"]][["thiosulfate->sulfate"]] = c("soxB", "soxC", "soxD")
    gene_check[["sulfur"]][["thiosulfate->tetrathionate-oxidation"]] = c("doxA", "doxD", "tsdA")
    gene_check[["sulfur"]][["tetrathionate->thiosulfate-reduction"]] = c("ttrA", "ttrB_dsrO", "ttrC")
    gene_check[["sulfur"]][["DMSO->DMS_and"]] = c("dmsA", "dmsB", "dmsC")
    gene_check[["sulfur"]][["DMSO->DMS_or"]] = c("dmsA", "dmsB", "dmsC")
    gene_check[["sulfur"]][["DMSP->DMS"]] = c("dddL")
    gene_check[["sulfur"]][["hdr"]] = c("hdrA1", "hdrB1", "hdrC1", "hdrA2", "hdrB2", "hdrC2")
    gene_check[["sulfur"]][["sulfate->sulfite-assimilatory"]] = c("cysI", "cysJ", "sir", "cysC", "cysH")
    gene_check[["sulfur"]][["sulfate->APS-assimilatory"]] = c("cysNC", "cysN", "cysD")

    gene_check[["carbonfixation"]][["rac"]] = c("cdhA","cdhB", "cooF", "cooS", "fdh.NADP+A", "fdh.NADP+B", "ftfL", "folD", "metF", "acsE", "acsB", "acsC", "acsD", "cooF", "cooS", "pta1", "pta2", "pta3", "ackA")

    gene_check[["phototrophy"]][["oxygenic photosystem I"]] = c("psaA", "psaB", "psaC", "psaD", "psaE", "psaF", "psaI", "psaJ", "psaK", "psaL", "psaM", "psaX")
    gene_check[["phototrophy"]][["oxygenic photosystem II"]] = c("psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbO", "psbP", "psbT", "psbY", "psbZ", "psb27", "psb28", "psbU", "psbV", "psbX", "psb28-2")
    gene_check[["phototrophy"]][["cytochrome b6f"]] = c("petA", "petB", "petC", "petD", "petG", "petL", "petM", "petN")
    gene_check[["phototrophy"]][["photosynthetic ETC"]] = c("petE", "petF", "petH", "petJ")
    gene_check[["phototrophy"]][["anoxygenic photosystem I RC"]] = c("pscA", "pscB", "pscC", "pscD", "fmoA", "csmA", "csmB", "csmC", "csmD", "csmE", "csmF", "csmH", "csmI", "csmJ", "csmX")
    gene_check[["phototrophy"]][["anoxygenic photosystem II RC"]] = c("pufL", "pufM", "puhA", "pufX", "pufC", "pufA", "pufB", "pucA", "pucB")
    gene_check[["phototrophy"]][["oxygenic chlorophyll a synthesis"]] = c("bchG")
    gene_check[["phototrophy"]][["bacteriochlorophyll a synthesis"]] =  c("bchH", "bchM", "bchE", "chlE", "bciB", "bchJ", "por", "chlL", "bchX", "bchF", "bchC", "bchG", "bchP")
    gene_check[["phototrophy"]][["bacteriochlorophyll b synthesis"]] =  c("bchH", "bchM", "bchE", "chlE", "bciB", "bchJ", "por", "chlL", "bchX", "bchF", "bchC", "bchG", "bchP")
    gene_check[["phototrophy"]][["bacteriochlorophyll c synthesis"]] =  c("bchH", "bchM", "bchE", "chlE", "bciB", "bchJ", "por", "chlL", "bciC", "bchV", "bchU", "bchK")
    gene_check[["phototrophy"]][["bacteriochlorophyll d synthesis"]] =  c("bchH", "bchM", "bchE", "chlE", "bciB", "bchJ", "por", "chlL", "bciC", "bchV", "bchK")
    gene_check[["phototrophy"]][["bacteriochlorophyll e synthesis"]] =  c("bchH", "bchM", "bchE", "chlE", "bciB", "bchJ", "por", "chlL", "bciC", "bchV", "bchU", "bciD")

    gene_check[["fumarate_respiration"]] = c("frdA", "frdB", "frdC", "frdD", "frd.NADHA", "frd.NADHB", "frd.NADHC", "frd.NADHD", "frd.NADHE")
    gene_check[["arsenic"]] = c("arrA", "arrB")
    gene_check[["selenium"]] = c("serA", "serB", "serC")
    gene_check[["iron"]] = c("mtrA", "mtrB", "mtrC", "mtrD", "mtrE", "mtrF", "omcA", "pioA", "pioB", "pioC", "foxE", "foxY", "foxZ", "cyc2")

    faprotax_check = list()
    faprotax_check[["dark_thiosulfate_oxidation"]] = "('soxB' | 'soxC' | 'soxD' | 'tsdA')"
    faprotax_check[["dark_sulfite_oxidation"]] = "('sorA' | 'sorB' | 'soeA' | 'soeB' | 'soeC' | 'aprA' | 'aprB')"
    faprotax_check[["chlorate_reducers"]] = "('serA' | 'serB' | 'serC')"
    faprotax_check[["sulfur_respiration"]] = "('sat' | 'aprA' | 'aprB' | 'sreA' | 'sreB' | 'sreC' | 'sqr')"
    faprotax_check[["sulfate_respiration"]] = "('sat')"

    result = list(rule_check = rule_check, gene_check = gene_check, faprotax_check = faprotax_check)
    return(result)
}
