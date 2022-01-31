#' Make sure binary traits are of type logical
#'
#' @param genomeset_results
#'
#' @import dplyr
#' @return genomeset_results
#' @export
convert_traitdatatype = function(genomeset_results, binarytype = "logical") {
  # to glimpse into data types
  # genomeset_results[["trait_matrixatgranularity3"]] %>% glimpse()
  granularities = c("1", "2", "3")
  for(g in 1:length(granularities)) {
    object = switch(granularities[g],
                    "1" = "trait_matrixatgranularity1",
                    "2" = "trait_matrixatgranularity2",
                    "3" = "trait_matrixatgranularity3"
    )
    binary_traits = traits_listbygranularity[[as.numeric(granularities[g])]] %>%
      dplyr::filter(`microtrait_trait-type` == "binary") %>%
      dplyr::pull(`microtrait_trait-name`) %>%
      as.character()
    count_traits = traits_listbygranularity[[as.numeric(granularities[g])]] %>%
      dplyr::filter(`microtrait_trait-type` == "count") %>%
      dplyr::pull(`microtrait_trait-name`) %>%
      as.character()
    if(binarytype == "logical") {
      genomeset_results[[as.numeric(granularities[g])]] = genomeset_results[[as.numeric(granularities[g])]] %>%
        dplyr::mutate_at(binary_traits, as.logical)
    }
    if(binarytype == "factor") {
      genomeset_results[[as.numeric(granularities[g])]] = genomeset_results[[as.numeric(granularities[g])]] %>%
        dplyr::mutate_at(binary_traits, as.factor)
    }
    if(binarytype == "special") {
      genomeset_results[[as.numeric(granularities[g])]] = genomeset_results[[as.numeric(granularities[g])]] %>%
        dplyr::mutate_at(binary_traits,
                         funs(case_when(. == 0 ~ 0.0,
                                        . == 1 ~ 15.13181))
                        )
    }
  }
  genomeset_results
}

#' Normalize count traits
#'
#' @param genomeset_results
#'
#' @import dplyr
#' @return genomeset_results
#' @export
trait.normalize <- function(genomeset_results,
                            normby = "Estimated Size") {
  # genomeset_results = genomeset_results_wmetadata
  granularities = c("1", "2", "3")
  for(g in 1:length(granularities)) {
    object = switch(granularities[g],
                    "1" = "trait_matrixatgranularity1",
                    "2" = "trait_matrixatgranularity2",
                    "3" = "trait_matrixatgranularity3"
                   )
    # will normalize count traits only
    count_traits = traits_listbygranularity[[as.numeric(granularities[g])]] %>%
      dplyr::filter(`microtrait_trait-type` == "count") %>%
      dplyr::pull(`microtrait_trait-name`) %>%
      as.character()
    if(normby != "zscore") {
      genomeset_results[[object]] = genomeset_results[[object]] %>%
        dplyr::filter(!is.na(get(normby))) %>% # shouldn't be common-there is one with no size information, taxon_id:2634166898
        dplyr::mutate_at(count_traits,
                         funs(./as.numeric(get(normby)) * 1E6))
    }
    if(normby == "zscore") {
      genomeset_results[[object]] = genomeset_results[[object]] %>%
        dplyr::mutate_at(count_traits,
                         funs((. - mean(.)) / sd(.)))
    }
    if(normby == "range") {
      genomeset_results[[object]] = genomeset_results[[object]] %>%
        dplyr::mutate_at(count_traits,
                         funs((. - min(.)) / (max(.) - min(.))))
    }
    # check calculations
    # i = 30
    # norm1 = genomeset_results_wmetadata[[object]] %>% dplyr::pull(count_traits[i]) /
    #   genomeset_results_wmetadata[[object]] %>% dplyr::pull(`Estimated Size`) *1E6
    # norm2 = genomeset_results[[object]] %>% dplyr::pull(count_traits[i])
    # table(norm1-norm2)
    # binary_traits = traits_listbygranularity[[as.numeric(granularities[g])]] %>%
    #   dplyr::filter(`microtrait_trait-type` == "binary") %>%
    #   dplyr::pull(`microtrait_trait-name`) %>%
    #   as.character()
    # table(genomeset_results_wmetadata[[object]] %>% dplyr::pull(binary_traits[i]) -
    #         genomeset_results[[object]] %>% dplyr::pull(binary_traits[i]))
  }
  genomeset_results
}

#' Convert continuous traits to binary
#'
#' @param feature_matrix
#' @import dplyr
#' @return results
#'
#' @export trait.continuous2binary
trait.continuous2binary <- function(feature_matrix) {
  feature_matrix_binary = feature_matrix %>%
    dplyr::mutate_at(vars(!starts_with(c("id", "mingentime", "ogt"))),
                     funs(case_when(. > 0 ~ 1,
                                    . <= 0 ~ 0))) # %>%
  # if prefer to create new columns
  #.funs = list(binary = ~case_when(. >= 1 ~ 1,TRUE ~ 0)))
  feature_matrix_binary
}

#' Fix traits
#'
#'
#'
# base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
# dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
# genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))
# genomeset_results_wmetadata = redefine.binarytraits(genomeset_results_wmetadata,trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate pathway")
# genomeset_results_wmetadata = redefine.binarytraits(genomeset_results_wmetadata,trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate/4-hydroxybutyrate pathway")
# genomeset_results_wmetadata = redefine.binarytraits(genomeset_results_wmetadata,trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:dicarboxylate-hydroxybutyrate pathway")
# genomeset_results_wmetadata = redefine.binarytraits(genomeset_results_wmetadata,trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive acetyl-CoA pathway")
# genomeset_results_wmetadata = redefine.binarytraits(genomeset_results_wmetadata,trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive citric acid pathway")
# genomeset_results_wmetadata = redefine.binarytraits(genomeset_results_wmetadata,trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:Calvin Cycle")
# genomeset_results_wmetadata = redefine.binarytraits(genomeset_results_wmetadata,trait = "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex III")
# genomeset_results_wmetadata = redefine.binarytraits(genomeset_results_wmetadata,trait = "Resource Use:Phototrophy:pigments:accessory:carotenoid:carotene:diaponeurosporene")
# saveRDS(genomeset_results_wmetadata, file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))
redefine.binarytraits <- function(genomeset_results_wmetadata,trait) {
  library(dplyr)
  rule_matrix = genomeset_results_wmetadata[["rule_matrix"]]
  hmm_matrix = genomeset_results_wmetadata[["hmm_matrix"]]
  trait_matrix3 = genomeset_results_wmetadata[["trait_matrixatgranularity3"]]

  #trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate pathway"
  #1. Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate pathway	-->	3hpb1_full
  if(trait == "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate pathway") {
    cat("Processing ", trait, "\n")
    ## reset
    trait_matrix3.1 = trait_matrix3 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate pathway`= 0)
    ## fix
    ids = rule_matrix %>% dplyr::filter(`malonyl-CoA->malonate semialdehyde`==TRUE | `acryloyl-CoA->propionyl-CoA`==TRUE) %>% dplyr::pull(`id`)
    trait_matrix3.1 = trait_matrix3.1 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate pathway` = ifelse(`id` %in% ids,1,0))
    genomeset_results_wmetadata[["trait_matrixatgranularity3"]] = trait_matrix3.1

    #trait_matrix3.1 %>%
    #  dplyr::select(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate pathway`) %>%
    #  table()
  }

  #trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate/4-hydroxybutyrate pathway"
  #2. Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate/4-hydroxybutyrate pathway	-->	3hp4hb_full
  if(trait == "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate/4-hydroxybutyrate pathway") {
    cat("Processing ", trait, "\n")
    ## reset
    trait_matrix3.1 = trait_matrix3 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate/4-hydroxybutyrate pathway`= 0)
    ## fix
    ids = rule_matrix %>% dplyr::filter((((`4-hydroxybutyryl-CoA->crotonoyl-CoA`==TRUE) &
                                            (`crotonoyl-CoA->3-hydroxybutyryl-CoA`==TRUE) &
                                            (`3-hydroxybutyryl-CoA->acetoacetyl-CoA`==TRUE) &
                                            (`acetoacetyl-CoA->acetyl-CoA`==TRUE)) |
                                           ((`succinyl-CoA->succinate semialdehyde`==TRUE) &
                                              (`succinate semialdehyde->4-hydroxybutyrate`==TRUE) &
                                              (`4-hydroxybutyryl-CoA->crotonoyl-CoA`==TRUE) &
                                              (`crotonoyl-CoA->3-hydroxybutyryl-CoA`==TRUE) &
                                              (`3-hydroxybutyryl-CoA->acetoacetyl-CoA`==TRUE) &
                                              (`acetoacetyl-CoA->acetyl-CoA`==TRUE) &
                                              (`acetyl-CoA+HCO3->malonyl-CoA`==TRUE) &
                                              (`acryloyl-CoA->propionyl-CoA`==TRUE)))) %>% dplyr::pull(`id`)
    trait_matrix3.1 = trait_matrix3.1 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate/4-hydroxybutyrate pathway` = ifelse(`id` %in% ids,1,0))
    genomeset_results_wmetadata[["trait_matrixatgranularity3"]] = trait_matrix3.1

    #trait_matrix3.1 %>%
    #  dplyr::filter(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:3-hydroxypropionate/4-hydroxybutyrate pathway` == 1) %>%
    #  dplyr::select(`NCBI_Phylum`) %>%
    #  table()
  }

  # trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:dicarboxylate-hydroxybutyrate pathway"
  #3. Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:dicarboxylate-hydroxybutyrate pathway	-->	dicarb-4hb_keyandCO2fixation
  if(trait == "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:dicarboxylate-hydroxybutyrate pathway") {
    cat("Processing ", trait, "\n")
    ## reset
    trait_matrix3.1 = trait_matrix3 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:dicarboxylate-hydroxybutyrate pathway`= 0)
    ## fix
    ids = rule_matrix %>% dplyr::filter(`dicarb-4hb_CO2fixation`==TRUE) %>% dplyr::pull(`id`)
    trait_matrix3.1 = trait_matrix3.1 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:dicarboxylate-hydroxybutyrate pathway` = ifelse(`id` %in% ids,1,0))
    genomeset_results_wmetadata[["trait_matrixatgranularity3"]] = trait_matrix3.1

    #trait_matrix3.1 %>%
    #  dplyr::filter(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:dicarboxylate-hydroxybutyrate pathway` == 1) %>%
    #  dplyr::select(`NCBI_Phylum`) %>%
    #  table()
  }

  # trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive acetyl-CoA pathway"
  #4. Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive acetyl-CoA pathway	-->	rac_full
  if(trait == "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive acetyl-CoA pathway") {
    cat("Processing ", trait, "\n")
    ## reset
    trait_matrix3.1 = trait_matrix3 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive acetyl-CoA pathway`= 0)
    ## fix
    ids = rule_matrix %>% dplyr::filter((`CO2->formate_or`==TRUE) | (`CO2->CO`==TRUE)) %>% dplyr::pull(`id`)
    trait_matrix3.1 = trait_matrix3.1 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive acetyl-CoA pathway` = ifelse(`id` %in% ids,1,0))
    genomeset_results_wmetadata[["trait_matrixatgranularity3"]] = trait_matrix3.1

    #trait_matrix3.1 %>%
    #  dplyr::filter(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive acetyl-CoA pathway` == 1) %>%
    #  dplyr::select(`NCBI_Phylum`) %>%
    #  table()
  }

  # trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive citric acid pathway"
  #5. Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive citric acid pathway	-->	rca_full
  if(trait == "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive citric acid pathway") {
    cat("Processing ", trait, "\n")
    ## reset
    trait_matrix3.1 = trait_matrix3 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive citric acid pathway`= 0)
    ## fix
    ids = rule_matrix %>% dplyr::filter(`citrate->oxaloacetate+acetyl-CoA_or`==TRUE) %>% dplyr::pull(`id`)
    trait_matrix3.1 = trait_matrix3.1 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive citric acid pathway` = ifelse(`id` %in% ids,1,0))
    genomeset_results_wmetadata[["trait_matrixatgranularity3"]] = trait_matrix3.1

    #trait_matrix3.1 %>%
    #  dplyr::filter(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:reductive citric acid pathway` == 1) %>%
    #  dplyr::select(`NCBI_Phylum`) %>%
    #  table()
  }

  #trait = "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:Calvin Cycle"
  # 6.
  if(trait == "Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:Calvin Cycle") {
    cat("Processing ", trait, "\n")
    ## reset
    trait_matrix3.1 = trait_matrix3 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:Calvin Cycle`= 0)
    ## fix
    ids = rule_matrix %>% dplyr::filter((`RuBP+CO2->3PG`==TRUE) | (`R5P->RuBP`==TRUE)) %>% dplyr::pull(`id`)
    trait_matrix3.1 = trait_matrix3.1 %>% dplyr::mutate(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:Calvin Cycle` = ifelse(`id` %in% ids,1,0))
    genomeset_results_wmetadata[["trait_matrixatgranularity3"]] = trait_matrix3.1

    #trait_matrix3.1 %>%
    #  dplyr::filter(`Resource Acquisition:Substrate assimilation:C1 compounds:CO2 fixation:Calvin Cycle` == 1) %>%
    #  dplyr::select(`NCBI_Phylum`) %>%
    #  table()
  }

  # 7. Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex III -->
  #trait="Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex III"
  if(trait == "Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex III") {
    cat("Processing ", trait, "\n")
    ## reset
    trait_matrix3.1 = trait_matrix3 %>% dplyr::mutate(`Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex III`= 0)
    ## fix
    ids = hmm_matrix %>%
      dplyr::filter(qcrA== 1 | fbcH==1 | petB_1==1 | petB_2==1 | qcrB==1) %>%
      dplyr::pull(`id`)
    trait_matrix3.1 = trait_matrix3.1 %>% dplyr::mutate(`Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex III` = ifelse(`id` %in% ids,1,0))
    genomeset_results_wmetadata[["trait_matrixatgranularity3"]] = trait_matrix3.1

    #trait_matrix3.1 %>%
    #  dplyr::filter(`Resource Use:Chemotrophy:chemoorganoheterotrophy:aerobic respiration:electron transport chain: ETC complex III` == 1) %>%
    #  dplyr::select(`NCBI_Phylum`) %>%
    #  table()
  }

  # 8. Resource Use:Phototrophy:pigments:accessory:carotenoid:carotene:diaponeurosporene -->
  #trait="Resource Use:Phototrophy:pigments:accessory:carotenoid:carotene:diaponeurosporene"
  if(trait == "Resource Use:Phototrophy:pigments:accessory:carotenoid:carotene:diaponeurosporene") {
    cat("Processing ", trait, "\n")
    ## reset
    trait_matrix3.1 = trait_matrix3 %>% dplyr::mutate(`Resource Use:Phototrophy:pigments:accessory:carotenoid:carotene:diaponeurosporene`= 0)
    ## fix
    ids = hmm_matrix %>%
      dplyr::filter(crtM== 1 | crtN==1) %>%
      dplyr::pull(`id`)
    trait_matrix3.1 = trait_matrix3.1 %>% dplyr::mutate(`Resource Use:Phototrophy:pigments:accessory:carotenoid:carotene:diaponeurosporene` = ifelse(`id` %in% ids,1,0))
    genomeset_results_wmetadata[["trait_matrixatgranularity3"]] = trait_matrix3.1

    #trait_matrix3.1 %>%
    #  dplyr::filter(`Resource Use:Phototrophy:pigments:accessory:carotenoid:carotene:diaponeurosporene` == 1) %>%
    #  dplyr::select(`NCBI_Phylum`) %>%
    #  table()
  }

  #trait_matrix3 %>%
  #  dplyr::inner_join(hmm_matrix, by = c("id" = "id")) %>%
  #  dplyr::filter((petE==1 & petF==1 & petH==1 & petJ==1)) %>%
  #  #dplyr::mutate(petEFHJ = paste0(petE,petF,petH,petJ)) %>%
  #  #dplyr::select(`petEFHJ`)%>%table
  #  dplyr::select(`Resource Use:Phototrophy:photosystem:oxygenic:photosynthetic electron transport`) = 1
  genomeset_results_wmetadata
}
