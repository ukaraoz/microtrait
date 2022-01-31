base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))

select = list.files("/opt/OGT_prediction-1.0.2_modified/genomes", pattern = ".fna") %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(id = sub(".fna", "", value)) %>%
  dplyr::select(`id`)

feature_matrix = genomeset_results_wmetadata[["trait_matrixatgranularity3"]]
for(r in 1:nrow(feature_matrix)) {
  #feature_matrix %>%
  #  dplyr::mutate(id_fa = paste0(`id`, ".fna")) %>%
  #  dplyr::select(c(`id_fa`, `id`)) %>%
  #  #dplyr::inner_join(select, by = c("id" = "id")) %>%
  #  write.table(row.names = F, col.names = F, quote = F, sep = "\t",
  #              #file = "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait/inst/extdata/microtrait_genomes_retrieved.txt")
  #              file = "/opt/OGT_prediction-1.0.2_modified/microtrait_genomes_retrieved.txt")
  outfile = paste0("/opt/OGT_prediction-1.0.2_modified/genomes/", feature_matrix[r,] %>% pull(`id`), "_OGTgenomeinfo.txt")
  cat("Writing ", outfile, "\n")
  feature_matrix[r,] %>%
    dplyr::mutate(id_fa = paste0(`id`, ".fna")) %>%
    dplyr::select(c(`id_fa`, `id`)) %>%
    write.table(row.names = F, col.names = F, quote = F, sep = "\t",
                file = outfile)
}

for(r in 1:nrow(feature_matrix)) {
  #feature_matrix %>%
  #  dplyr::select(c(`id`, `NCBI_Superkingdom`, `NCBI_Phylum`, `NCBI_Class`, `NCBI_Order`, `NCBI_Family`)) %>%
  #  dplyr::rename(species = id,
  #                superkingdom = NCBI_Superkingdom,
  #                phylum = NCBI_Phylum,
  #                class = NCBI_Class,
  #                order = NCBI_Order,
  #                family = NCBI_Family) %>%
  #  #dplyr::inner_join(select, by = c("species" = "id")) %>%
  #  write.table(row.names = F, col.names = T, quote = F, sep = "\t",
  #              #file = "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait/inst/extdata/microtrait_species_taxonomic.txt")
  #              file = "/opt/OGT_prediction-1.0.2_modified/microtrait_species_taxonomic.txt")
  outfile = paste0("/opt/OGT_prediction-1.0.2_modified/genomes/", feature_matrix[r,] %>% pull(`id`), "_OGTspeciestaxonomicinfo.txt")
  cat(r, ": Writing ", outfile, "\n")
  feature_matrix[r,] %>%
    dplyr::select(c(`id`, `NCBI_Superkingdom`, `NCBI_Phylum`, `NCBI_Class`, `NCBI_Order`, `NCBI_Family`)) %>%
    dplyr::rename(species = id,
                  superkingdom = NCBI_Superkingdom,
                  phylum = NCBI_Phylum,
                  class = NCBI_Class,
                  order = NCBI_Order,
                  family = NCBI_Family) %>%
    #dplyr::inner_join(select, by = c("species" = "id")) %>%
    write.table(row.names = F, col.names = T, quote = F, sep = "\t",
                #file = "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait/inst/extdata/microtrait_species_taxonomic.txt")
                file = outfile)
}
