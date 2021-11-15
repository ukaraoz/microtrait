library(microtrait)
library(dplyr)
library(ggplot2)
library(ape)

# defaults
local = FALSE
subsample = TRUE
subsample.frac = 0.5
isolates = TRUE

if(local) {
  base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
  outdir = file.path(base, "traitconservation/consentrait_out")
  consentraitbin = "/Users/ukaraoz/Work/microtrait/code/microtrait/scripts/consentrait.R"
} else {
  base = "~/cscratch/alltarballs"
  outdir = file.path(base, "traitconservation/consentrait_out1")
  consentraitbin = "/global/homes/u/ukaraoz/bin/consentrait.R"
}

dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))
tree_file = file.path(base, "raxml_rpsC.rooted.tree")
tree = read.tree(tree_file)
trait_matrix = genomeset_results_wmetadata[["trait_matrixatgranularity3"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Aquatic", "Terrestrial", "Plants"))
#allids =trait_matrix %>% pull(`id`)
traits = c(traits_listbygranularity[[3]] %>%
             #dplyr::filter(`microtrait_trait-type` %in% c("count","count_by_substrate")) %>%
             dplyr::select(`microtrait_trait-name`) %>% pull() %>% as.character())

# binarize
trait_matrix = trait_matrix %>% dplyr::select(c(`id`, traits)) %>%
  dplyr::mutate_at(vars(!starts_with("id")),
                   funs(case_when(. >= 1 ~ 1,TRUE ~ 0)))

if(subsample == TRUE) {
  ngenomes = floor(nrow(trait_matrix)*subsample.frac)
} else {
  ngenomes = nrow(trait_matrix)
}

if(isolates == TRUE) {
  treeroot = "2757320680_2758412355_Archaea_Candidatus"
  trait_matrix.sample = dplyr::bind_rows(trait_matrix %>% sample_n(ngenomes),
                                         trait_matrix %>% dplyr::filter(`id` == "2757320680")
  )
} else {
  treeroot = "Archaea_QMVG01000001.1_100__GCA_003660975.1_ASM366097v1"
  trait_matrix.sample = dplyr::bind_rows(trait_matrix %>% sample_n(ngenomes),
                                         trait_matrix %>% dplyr::filter(`id` == "Archaea_QMVG01000001.1_100__GCA_003660975.1_ASM366097v1")
  )
}

# add matching tree tips
trait_matrix.sample = tree$tip.label %>% tibble::as_tibble() %>%
  dplyr::rename(treetip = value) %>%
  dplyr::mutate(taxonid = stringr::str_replace(`treetip`, "_.*", "")) %>%
  dplyr::right_join(trait_matrix.sample, by = c("taxonid" = "id")) %>%
  dplyr::rename(`id` = `taxonid`) %>%
  dplyr::filter(!is.na(treetip))

# prune the tree to the subsample
tree.sample = keep.tip(tree, c(treeroot, trait_matrix.sample %>% pull(`treetip`)))
for(i in 1:30) {
  cat(traits[i], "\n")

  status = try(run.consentrait.singletrait(traits[i],
                                  traitmatrix = trait_matrix.sample,
                                  tree = tree.sample,
                                  treeroot = treeroot,
                                  nbootstrap = 100,
                                  nproc = 40,
                                  perc.share.cutoff = 90,
                                  outdir = outdir,
                                  consentraitbin = consentraitbin,
                                  run = TRUE),
              silent = TRUE)
  if(class(status) == "try-error") {
    status = data.frame(command = traits[i], time = "error")
  }
}
