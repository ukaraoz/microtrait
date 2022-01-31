library(dplyr)
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
