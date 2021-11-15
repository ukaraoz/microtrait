library(microtrait)
library(dplyr)
library(ggplot2)
library(nlme)
library(ape)
#base = "/global/homes/u/ukaraoz/cscratch/alltarballs"
base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"

dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))

trait_matrix = genomeset_results_wmetadata[["trait_matrixatgranularity3"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Aquatic", "Terrestrial", "Plants"))

traits = c("mingentime", "OGT", traits_listbygranularity[[3]] %>%
  #dplyr::filter(`microtrait_trait-type` %in% c("count","count_by_substrate")) %>%
  dplyr::select(`microtrait_trait-name`) %>% pull() %>% as.character())

models = list()
for(i in 1:length(traits)) {
  cat(i, "\n")

  models[[i]] = build.taxa.model(trait_matrix, trait = traits[i], species = TRUE, optimizer = "nlm")
}
saveRDS(models, file = file.path(base, "traitconservation", "nlme.models.optimizer-nlm.species.rds"))

# ncores = 20
# models = parallel::mclapply(1:length(traits),
#                             function(i) {
#                               returnList = build.taxa.model(trait_matrix, trait = traits[i])
#                               returnList
#                             },
#                             mc.cores = ncores)


