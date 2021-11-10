library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(usethis)

base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))[["trait_matrixatgranularity3"]]

faprotax_taxa2ncbi_taxa = readr::read_delim("./data-raw/faprotax_taxa2ncbi_taxa.txt", delim = "\t")
message("Number of faprotax taxa: ", faprotax_taxa2ncbi_taxa %>% dplyr::select(`faprotax_taxa`) %>% n_distinct(), "\n")
#clade   class  family   genus   order  phylum species
#1       3      42     927      21       4    3241
#message("Number of faprotax species: ", faprotax_taxa2ncbi_taxa %>% dplyr::filter(`rank` == "species") %>% dplyr::select(`faprotax_taxa`) %>% n_distinct(), "\n")

# 7,831
faprotax = readr::read_delim("./data-raw/faprotax.txt", quote = "", delim = "\t") %>%
  dplyr::filter(`relationtype` == "taxa") %>%
  dplyr::left_join(faprotax_taxa2ncbi_taxa, by = c("taxa" = "faprotax_taxa"))

# find matching genomes for faprotax for each taxonomic rank
## species: 2,288
faprotax_species2genomes = faprotax_taxa2ncbi_taxa %>% dplyr::filter(`rank` == "species") %>% dplyr::inner_join(genomeset_results_wmetadata, by = c("ncbi_taxa" = "NCBI_Species")) %>% dplyr::select(c("faprotax_taxa", "ncbi_taxa", "rank", "id"))
  #dplyr::filter(`faprotax_taxa` == "Methanosarcina mazei") #dplyr::group_by(`faprotax_taxa`) %>% dplyr::tally() %>% dplyr::arrange(n)
## genus: 6,065
faprotax_genus2genomes = faprotax_taxa2ncbi_taxa %>% dplyr::filter(`rank` == "genus") %>% dplyr::inner_join(genomeset_results_wmetadata, by = c("ncbi_taxa" = "NCBI_Genus")) %>% dplyr::select(c("faprotax_taxa", "ncbi_taxa", "rank", "id"))
## family:  1,253
faprotax_family2genomes = faprotax_taxa2ncbi_taxa %>%  dplyr::filter(`rank` == "family") %>% dplyr::inner_join(genomeset_results_wmetadata, by = c("ncbi_taxa" = "NCBI_Family")) %>% dplyr::select(c("faprotax_taxa", "ncbi_taxa", "rank", "id"))
## order:   1,055
faprotax_order2genomes = faprotax_taxa2ncbi_taxa %>%  dplyr::filter(`rank` == "order") %>% dplyr::inner_join(genomeset_results_wmetadata, by = c("ncbi_taxa" = "NCBI_Order")) %>% dplyr::select(c("faprotax_taxa", "ncbi_taxa", "rank", "id"))
## class:   35
faprotax_class2genomes = faprotax_taxa2ncbi_taxa %>%  dplyr::filter(`rank` == "class") %>% dplyr::inner_join(genomeset_results_wmetadata, by = c("ncbi_taxa" = "NCBI_Class")) %>% dplyr::select(c("faprotax_taxa", "ncbi_taxa", "rank", "id"))
## phylum:  27
faprotax_phylum2genomes = faprotax_taxa2ncbi_taxa %>%  dplyr::filter(`rank` == "phylum") %>% dplyr::inner_join(genomeset_results_wmetadata, by = c("ncbi_taxa" = "NCBI_Phylum")) %>% dplyr::select(c("faprotax_taxa", "ncbi_taxa", "rank", "id"))
# faprotax_taxa2ncbi_taxa %>%  dplyr::filter(`rank` == "nonspecies")
faprotax2genomes = dplyr::bind_rows(faprotax_phylum2genomes, faprotax_class2genomes, faprotax_order2genomes, faprotax_family2genomes, faprotax_genus2genomes, faprotax_species2genomes)

usethis::use_data(faprotax, overwrite = TRUE)
usethis::use_data(faprotax2genomes, overwrite = TRUE)
#faprotax_species = faprotax %>% dplyr::filter(`rank` == "species")  #2706 species
#
#genomeset_results_wmetadata %>% dplyr::select(c("id", "NCBI_Species")
#
#faprotax_species_unmatchedgenomes = faprotax_species %>%
#  dplyr::anti_join(genomeset_results_wmetadata, by = c("ncbi_taxa" = "NCBI_Species"))


