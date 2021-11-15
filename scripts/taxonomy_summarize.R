library(tidyr)
library(dplyr)

# summarize IMG genomes taxonomy

base = "/Users/ukaraoz/Work/microtrait/code"
ranks = c("NCBI Phylum", "NCBI Class","NCBI Order", "NCBI Family", "NCBI Genus", "NCBI Species")
img_genomes = readRDS(file.path(base, "inst/extdata/", "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID.metadata.rds"))
img_genomes = img_genomes %>%
  dplyr::filter((`Ecosystem` == "Environmental"|`Ecosystem` == "Host-associated")&
                (`Ecosystem Category`%in%c("Aquatic", "Terrestrial", "Plants"))) %>%
  tidyr::unite(Ecosystem,"Ecosystem Category", "Ecosystem Type", sep = ":", remove = FALSE) %>%
  dplyr::select(c("IMG Taxon ID", "NCBI Superkingdom", "NCBI Phylum", "NCBI Class","NCBI Order",
                  "NCBI Family", "NCBI Genus", "NCBI Species", "Ecosystem Category", "Ecosystem Type", "Ecosystem"))

taxonids=img_genomes %>% dplyr::pull(`IMG Taxon ID`)

for(i in 1:length(ranks)) {
  result = summarize_taxonomy(taxonids, ranks[i])
  write.table(result %>% as.data.frame(),
              file = file.path(base, "inst/extdata/", paste0(sub("NCBI ", "img_genomes_", ranks[i]), ".stats.xls")),
              quote = F, row.names = F, col.names = T, sep = "\t")
}

summarize_taxonomy <- function(taxonids, rank) {
  if(rank == "NCBI Phylum") {
    temp = img_genomes %>%
      dplyr::filter(`IMG Taxon ID` %in% taxonids) %>%
      dplyr::mutate(rank = get(rank)) %>%
      tidyr::unite(taxa, "NCBI Superkingdom", "rank", sep = ":") %>%
      dplyr::select(c("IMG Taxon ID", taxa, "Ecosystem Category")) %>%
      #tidyr::pivot_longer(!`IMG Taxon ID`, names_to = "column", values_to = "value") %>%
      #dplyr::select(-`IMG Taxon ID`) %>%
      #dplyr::count(column, value,.drop=FALSE) %>%
      dplyr::count(`taxa`, `Ecosystem Category`) %>%
      dplyr::mutate(percent = n/sum(n)*100)
      #dplyr::mutate(rank = as.factor(rank),
      #              `Ecosystem Category` = as.factor(`Ecosystem Category`))
  }
  if(rank == "NCBI Class") {
    temp = img_genomes %>%
      dplyr::filter(`IMG Taxon ID` %in% taxonids) %>%
      dplyr::mutate(rank = get(rank)) %>%
      tidyr::unite(taxa, "NCBI Superkingdom", "NCBI Phylum", "rank", sep = ":") %>%
      dplyr::select(c("IMG Taxon ID", taxa, "Ecosystem Category")) %>%
      #tidyr::pivot_longer(!`IMG Taxon ID`, names_to = "column", values_to = "value") %>%
      #dplyr::select(-`IMG Taxon ID`) %>%
      #dplyr::count(column, value,.drop=FALSE) %>%
      dplyr::count(`taxa`, `Ecosystem Category`) %>%
      dplyr::mutate(percent = n/sum(n)*100)
      #dplyr::mutate(rank = as.factor(rank),
      #              `Ecosystem Category` = as.factor(`Ecosystem Category`))
  }
  if(rank == "NCBI Order") {
    temp = img_genomes %>%
      dplyr::filter(`IMG Taxon ID` %in% taxonids) %>%
      dplyr::mutate(rank = get(rank)) %>%
      tidyr::unite(taxa, "NCBI Superkingdom", "NCBI Phylum", , "NCBI Class", "rank", sep = ":") %>%
      dplyr::select(c("IMG Taxon ID", taxa, "Ecosystem Category")) %>%
      #tidyr::pivot_longer(!`IMG Taxon ID`, names_to = "column", values_to = "value") %>%
      #dplyr::select(-`IMG Taxon ID`) %>%
      #dplyr::count(column, value,.drop=FALSE) %>%
      dplyr::count(`taxa`, `Ecosystem Category`) %>%
      dplyr::mutate(percent = n/sum(n)*100)
      #dplyr::mutate(rank = as.factor(rank),
      #              `Ecosystem Category` = as.factor(`Ecosystem Category`))
  }
  if(rank == "NCBI Family") {
    temp = img_genomes %>%
      dplyr::filter(`IMG Taxon ID` %in% taxonids) %>%
      dplyr::mutate(rank = get(rank)) %>%
      tidyr::unite(taxa, "NCBI Superkingdom", "NCBI Phylum", "NCBI Class", "NCBI Order", "rank", sep = ":") %>%
      dplyr::select(c("IMG Taxon ID", taxa, "Ecosystem Category")) %>%
      #tidyr::pivot_longer(!`IMG Taxon ID`, names_to = "column", values_to = "value") %>%
      #dplyr::select(-`IMG Taxon ID`) %>%
      #dplyr::count(column, value,.drop=FALSE) %>%
      dplyr::count(`taxa`, `Ecosystem Category`) %>%
      dplyr::mutate(percent = n/sum(n)*100)
      #dplyr::mutate(rank = as.factor(rank),
      #              `Ecosystem Category` = as.factor(`Ecosystem Category`))
  }
  if(rank == "NCBI Genus") {
    temp = img_genomes %>%
      dplyr::filter(`IMG Taxon ID` %in% taxonids) %>%
      dplyr::mutate(rank = get(rank)) %>%
      tidyr::unite(taxa, "NCBI Superkingdom", "NCBI Phylum", "NCBI Class", "NCBI Order", "NCBI Family", "rank", sep = ":") %>%
      dplyr::select(c("IMG Taxon ID", taxa, "Ecosystem Category")) %>%
      #tidyr::pivot_longer(!`IMG Taxon ID`, names_to = "column", values_to = "value") %>%
      #dplyr::select(-`IMG Taxon ID`) %>%
      #dplyr::count(column, value,.drop=FALSE) %>%
      dplyr::count(`taxa`, `Ecosystem Category`) %>%
      dplyr::mutate(percent = n/sum(n)*100)
      #dplyr::mutate(rank = as.factor(rank),
      #              `Ecosystem Category` = as.factor(`Ecosystem Category`))
  }
  if(rank == "NCBI Species") {
    temp = img_genomes %>%
      dplyr::filter(`IMG Taxon ID` %in% taxonids) %>%
      dplyr::mutate(rank = get(rank)) %>%
      tidyr::unite(taxa, "NCBI Superkingdom", "NCBI Phylum", "NCBI Class", "NCBI Order", "NCBI Family", "NCBI Genus", "rank", sep = ":") %>%
      dplyr::select(c("IMG Taxon ID", taxa, "Ecosystem Category")) %>%
      #tidyr::pivot_longer(!`IMG Taxon ID`, names_to = "column", values_to = "value") %>%
      #dplyr::select(-`IMG Taxon ID`) %>%
      #dplyr::count(column, value,.drop=FALSE) %>%
      dplyr::count(`taxa`, `Ecosystem Category`) %>%
      dplyr::mutate(percent = n/sum(n)*100)
      #dplyr::mutate(rank = as.factor(rank),
      #              `Ecosystem Category` = as.factor(`Ecosystem Category`))
  }
  percent = temp %>%
    dplyr::select(-n) %>%
    tidyr::pivot_wider(names_from = "Ecosystem Category",
                       values_from = "percent",
                       values_fill = 0) %>%
    dplyr::rename(`Aquatic_percent` = `Aquatic`,
                  `Terrestrial_percent` = `Terrestrial`,
                  `Plant-associated_percent` = `Plants`)
  count = temp %>%
    dplyr::select(-percent) %>%
    tidyr::pivot_wider(names_from = "Ecosystem Category",
                       values_from = "n",
                       values_fill = 0) %>%
    dplyr::rename(`Aquatic_count` = `Aquatic`,
                  `Terrestrial_count` = `Terrestrial`,
                  `Plant-associated_count` = `Plants`)
  result = count %>% dplyr::left_join(percent, by = c("taxa" = "taxa"))
  return(result)
}
