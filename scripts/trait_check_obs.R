library(ggplot2)

base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))
objects = c("trait_matrixatgranularity3", "hmm_matrix", "rule_matrix")
for(i in 1:length(objects)) {
  feature_matrix = genomeset_results_wmetadata[[objects[i]]] %>%
    dplyr::filter(`Ecosystem_Category` %in% c("Aquatic", "Terrestrial", "Plants"))
  write.table(feature_matrix,
              file = file.path(base, paste0(dataset, ".", objects[i], ".xls")),
              row.names = F, col.names = T, sep = "\t", quote = F)

  feature_matrix_prevalence = compute.prevalence(feature_matrix, objects[i])
  write.table(feature_matrix_prevalence,
            file = file.path(base, paste0(dataset, ".", objects[i], "_prevalence.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)
}

i=3
feature_matrix = genomeset_results_wmetadata[[objects[i]]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Aquatic", "Terrestrial", "Plants"))
select_cols = c("id", "NCBI_Superkingdom", "NCBI_Phylum", "NCBI_Class", "NCBI_Order", "NCBI_Family", "NCBI_Genus")

# dicarb-4hb_keyandCO2fixation	NEW dicarb-4hb_CO2fixation


#select_rules_3hpb1_full = c("acetyl-CoA+HCO3->malonyl-CoA",
#"malonyl-CoA->malonate semialdehyde",
#"malonate semialdehyde->3-hydroxypropanoate",
#"3-hydroxypropanoate->3-hydroxypropanoyl-CoA",
#"3-hydroxypropanoyl-CoA->acryloyl-CoA",
#"acryloyl-CoA->propionyl-CoA",
#"propionyl-CoA+HCO3->S-methylmalonyl-CoA",
#"S-methylmalonyl-CoA<->R-methylmalonyl-CoA",
#"R-methylmalonyl-CoA<->succinyl-CoA",
#"malate->malyl-CoA_or",
#"succinate->fumarate_key",
#"fumarate->malate_full",
#"malyl-CoA->glyoxylate+acetyl-CoA",
#"propionyl-CoA->acetyl-CoA_glyoxylateassim",
#"propionyl-CoA+glyoxylate->beta-methylmalyl-CoA",
#"beta-methylmalyl-CoA->3-methylfumaroyl-CoA",
#"3-methylfumaroyl-CoA->3-methylfumaryl-CoA",
#"3-methylfumaryl-CoA->3S-citramalyl-CoA",
#"3S-citramalyl-CoA->acetyl-CoA+pyruvate")
# NEW 		malonyl-CoA->malonate semialdehyde OR acryloyl-CoA->propionyl-CoA
#feature_matrix_select = feature_matrix %>% select(c(select_cols, select_rules_3hpb1_full)) %>%
#  mutate_if(is.logical, as.integer)
#write.table(feature_matrix_select,
#            file = file.path(base, paste0(dataset, ".", objects[i], "_3hpb1.xls")),
#            row.names = F, col.names = T, sep = "\t", quote = F)

#select_rules_3hp4hb_full = c("propionyl-CoA+HCO3->S-methylmalonyl-CoA",
#                 "S-methylmalonyl-CoA<->R-methylmalonyl-CoA",
#                 "R-methylmalonyl-CoA<->succinyl-CoA",
#                 "succinyl-CoA->succinate semialdehyde",
#                 "succinate semialdehyde->4-hydroxybutyrate",
#                 "4-hydroxybutyrate->4-hydroxybutyryl-CoA",
#                 "4-hydroxybutyryl-CoA->crotonoyl-CoA",
#                 "crotonoyl-CoA->3-hydroxybutyryl-CoA",
#                 "3-hydroxybutyryl-CoA->acetoacetyl-CoA",
#                 "acetoacetyl-CoA->acetyl-CoA",
#                 "propionyl-CoA->acetyl-CoA-3hp4hb",
#                 "acetyl-CoA->propionyl-CoA",
#                 "acetyl-CoA+HCO3->malonyl-CoA",
#                 "malonyl-CoA->malonate semialdehyde",
#                 "malonate semialdehyde->3-hydroxypropanoate",
#                 "3-hydroxypropanoate->3-hydroxypropanoyl-CoA",
#                 "3-hydroxypropanoyl-CoA->acryloyl-CoA",
#                 "acryloyl-CoA->propionyl-CoA")
# 3hp4hb_full		NEW 4-hydroxybutyryl-CoA->crotonoyl-CoA & crotonoyl-CoA->3-hydroxybutyryl-CoA & 3-hydroxybutyryl-CoA->acetoacetyl-CoA &acetoacetyl-CoA->acetyl-CoA
# OR
# succinyl-CoA->succinate semialdehyde & succinate semialdehyde->4-hydroxybutyrate & 4-hydroxybutyryl-CoA->crotonoyl-CoA & crotonoyl-CoA->3-hydroxybutyryl-CoA & 3-hydroxybutyryl-CoA->acetoacetyl-CoA & acetoacetyl-CoA->acetyl-CoA &
#   acetyl-CoA+HCO3->malonyl-CoA & acryloyl-CoA->propionyl-CoA
#feature_matrix_select = feature_matrix %>% select(c(select_cols, select_rules_3hp4hb_full)) %>%
#  mutate_if(is.logical, as.integer)
#write.table(feature_matrix_select,
#                       file = file.path(base, paste0(dataset, ".", objects[i], "_3hp4hb.xls")),
#                       row.names = F, col.names = T, sep = "\t", quote = F)

#rac_full	NEW: ('CO2->formate_or') | ('CO2->CO')

#rca_full	NEW: citrate->oxaloacetate+acetyl-CoA_or

hmm_matrix = genomeset_results_wmetadata[["hmm_matrix"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Aquatic", "Terrestrial", "Plants"))
rules_matrix = genomeset_results_wmetadata[["rule_matrix"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Aquatic", "Terrestrial", "Plants"))

select_rules_cc = c("RuBP+CO2->3PG", "PGAP<->3PG","G3P<->PGAP", "DHAP<->G3P",
"G3P+DHAP->FBP","FBP->F6P","SDP->S7P","X5P+E4P<->F6P+G3P", "X5P->R5P","R5P->RuBP")
rules_matrix_select = feature_matrix %>% select(c(select_cols, select_rules_cc)) %>%
  mutate_if(is.logical, as.integer)
write.table(rules_matrix_select,
            file = file.path(base, paste0(dataset, ".", objects[i], "_cc.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

# ~ok
carboxysome_genes = c("eutK","eutL","eutM","eutN","eutS",
                      "ccmK","ccmL","ccmM","ccmN","ccmO")
hmm_matrix_select = hmm_matrix %>% select(c(select_cols, carboxysome_genes)) %>%
  mutate_if(is.logical, as.integer)
write.table(hmm_matrix_select,
            file = file.path(base, paste0(dataset, ".", objects[i], "_carboxysome_genes.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

# ok
anaerobicammoniaoxidation_genes = c("hzsA", "hzsB", "hzsC","hdh")
hmm_matrix_select = hmm_matrix %>% select(c(select_cols, anaerobicammoniaoxidation_genes)) %>%
  mutate_if(is.logical, as.integer)
write.table(hmm_matrix_select,
            file = file.path(base, paste0(dataset, ".", objects[i], "_anaerobicammoniaoxidation_genes.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

# CHANGE TO ('qcrA' | 'fbcH' | 'petB_2' | 'petB_1' | 'qcrB')
ETC_complexIII_genes = c("qcrA", "fbcH", "petB_2", "petB_1", "qcrB")
hmm_matrix_select = hmm_matrix %>% select(c(select_cols, ETC_complexIII_genes)) %>%
  mutate_if(is.logical, as.integer)
write.table(hmm_matrix_select,
            file = file.path(base, paste0(dataset, ".", objects[i], "_ETC_complexIII_genes.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

# leave it as is
ETC_complexIV_rules = c("cytochromeoubiquinol_and", "cytochromec_and",
                        "cytochromeaa3600menaquinol_and", "cytochromeccbb3oxidase_and")
rules_matrix_select = rules_matrix %>% select(c(select_cols, ETC_complexIV_rules)) %>%
  mutate_if(is.logical, as.integer)
write.table(rules_matrix_select,
            file = file.path(base, paste0(dataset, ".", objects[i], "_ETC_complexIV_genes.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

anoxygenicphotosystemIRC_genes = c("pscA", "pscB", "pscC", "pscD")
hmm_matrix_select = hmm_matrix %>% select(c(select_cols, anoxygenicphotosystemIRC_genes)) %>%
  mutate_if(is.logical, as.integer)
write.table(hmm_matrix_select,
            file = file.path(base, paste0(dataset, ".", objects[i], "_anoxygenicphotosystemIRC_genes.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

anoxygenicphotosystemIIRC_genes = c("pufL", "pufM", "puhA", "pufX", "pufC")
hmm_matrix_select = hmm_matrix %>% select(c(select_cols, anoxygenicphotosystemIIRC_genes)) %>%
  mutate_if(is.logical, as.integer)
write.table(hmm_matrix_select,
            file = file.path(base, paste0(dataset, ".", objects[i], "_anoxygenicphotosystemIIRC_genes.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

diaponeurosporenebiosynthesis_genes = c("crtM", "crtN")
hmm_matrix_select = hmm_matrix %>% select(c(select_cols, diaponeurosporenebiosynthesis_genes)) %>%
  mutate_if(is.logical, as.integer)
write.table(hmm_matrix_select,
            file = file.path(base, paste0(dataset, ".", objects[i], "_diaponeurosporenebiosynthesis_genes.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

#############################
# keyword = "denitrification"
# rulename2boolean = traits_listbygranularity[[2]] %>%
#   dplyr::filter(stringr::str_detect(`microtrait_trait-name`, keyword)) %>%
#   dplyr::left_join(rule2trait, by = c("microtrait_trait-name" = "microtrait_trait-name2")) %>%
#   dplyr::select(c("microtrait_rule-name", "microtrait_rule-boolean")) %>% as.data.frame()
# 
# ruleparts = data.frame()
# for(r in 1:nrow(rulename2boolean)) {
#   rulepartstemp= rule2parts(a[r,2])
#   temp = data.frame(`trait-keyword` = rep(keyword, length(rulepartstemp)),
#                     `microtrait_rule-name` = rep(as.character(rulename2boolean[r,1]), length(rulepartstemp)),
#                     `microtrait_rule-name-part` = rulepartstemp,
#                     check.names = FALSE)
#   ruleparts = rbind(ruleparts, temp)
# }

rule_check = list()
rule_check[["denitrification"]][["steps"]][["full"]] = c("nitrate->nitrite_dissimilatory_full",
                                                         "nitrite->nitricoxide_full",
                                                         "nitricoxide->nitrousoxide_full",
                                                         "nitrousoxide->dinitrogen")
rule_check[["denitrification"]][["steps"]][["partial"]] = c("nitrate->nitrite_dissimilatory_partial",
                                                            "nitrite->nitricoxide_partial",
                                                            "nitricoxide->nitrousoxide_partial",
                                                            "nitrousoxide->dinitrogen")
rule_check[["denitrification"]][["endproduct"]] = c("nitrate->nitrite_dissimilatory_full",
                                                    "nitrate->nitricoxide_complete",
                                                    "nitrate->nitrousoxide_complete",
                                                    "nitrate->dinitrogen_complete")
rule_check[["nitrogen fixation"]][["full"]] = "dinitrogen->ammonia_full"
rule_check[["nitrogen fixation"]][["partial"]] = "dinitrogen->ammonia_partial"
rule_check[["nitrogen fixation"]][["key"]] = "dinitrogen->ammonia_key"
rule_check[["ammonia oxidation"]][["full"]] = c("ammonia->hydroxylamine_full", "hydroxylamine->nitrite")
rule_check[["ammonia oxidation"]][["partial"]] = c("ammonia->hydroxylamine_partial", "hydroxylamine->nitrite")
rule_check[["nitrite oxidation"]][["full"]] = c("nitrite oxidation_full")
rule_check[["nitrite oxidation"]][["partial"]] = c("nitrite oxidation_partial")
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
rule_check[["ETC"]][["ETC_complexI"]][["full"]] = c("nuo_and")
rule_check[["ETC"]][["ETC_complexII"]][["full"]] = c("succinate->fumarate_full")
rule_check[["ETC"]][["ETC_complexIII"]][["full"]] = c("cytochromecreductase_or")
rule_check[["ETC"]][["ETC_complexIV"]][["full"]] = c("cytochromeoubiquinol_and", 
                                                     "cytochromec_and", 
                                                     "cytochromeaa3600menaquinol_and", 
                                                     "cytochromeccbb3oxidase_and")




  ETC_rules = c("nuo_and", "succinate->fumarate_key", "cytochromecreductase_or", "cytochromeoubiquinol_and", "cytochromec_and", "cytochromeaa3600menaquinol_and", "cytochromeccbb3oxidase_and")


gene_check = list()
gene_check[["denitrification"]] = c('narG', 'narH', 'narI_dsrM', 'napA', 'napB', 'nirS', 'nirK', 'norB', 'norC', 'norV', 'norW', 'nosZ')
gene_check[["nitrogen fixation"]] = c('nifD', 'nifK', 'nifH', 'anfG', 'vnfD', 'vnfG', 'vnfK', 'vnfH')
gene_check[["ammonia oxidation"]] = c('pmoA-amoA', 'pmoB-amoB', 'pmoC-amoC', 'hao')
gene_check[["nitrite oxidation"]] = c('narG' & 'narH')
gene_check[["methanogenesis"]] = c('tmtrA', 'tmtrB', 'tmtrC', 'tmtrD', 'tmtrE', 'tmtrF', 'tmtrG', 'tmtrH', 'mcrA', 'mcrB', 'mcrG', 'mcr2', 'mcrC', 'mcrD')
gene_check[["methanogenesis"]][["acetoclastic"]] = c('ackA', 'pta1', 'pta2', 'cdhC', 'acsC', 'acsD', gene_check[["methanogenesis"]])
gene_check[["methanogenesis"]][["hydrogenotrophic"]] = c('fwdA', 'fwdB', 'fwdC', 'fwdD', 'fwdF', 'fwdH', 'fwdG', 'fwdE', 'ftr', 'mtdch', 'hmd', 'mer', gene_check[["methanogenesis"]])
gene_check[["methanogenesis"]][["methylotrophic"]] = c('mmtsA', 'mtaB', 'mtaC', 'torA', 'torZ', 'torC', 'torY', 'mtbA', 'mtmB', 'mtmC', 'mtbA', 'mtbB', 'mtbC', 'mtbA', 'mttB', 'mttC', gene_check[["methanogenesis"]])

