library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ape)

base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
tree.file = file.path(base, "raxml_rpsC/RAxML_bestTree.result")
dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genome_metadata = readRDS(file.path(base, paste0(dataset, ".metadata.rds")))
bioproject2magid2domain = read.table(file.path(base, "bioproject2magid2domain.txt"), header = T, sep = "\t") %>% tibble::as_tibble()
treetip2rpsC = readr::read_tsv(file.path(base, "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID_treetip2rpsC.txt"),
                               col_types = cols(treetip = col_character(),
                                                geneid = col_character(),
                                                taxonid = col_character()
                                                )
                               ) %>% dplyr::select(c("geneid", "taxonid"))
tree = read.tree(tree.file)
# add domain names to tree tips
tree.tips.rename = tree$tip.label %>% tibble::as_tibble() %>%
  dplyr::rename(treetip = value) %>%
  tidyr::separate(`treetip`, c("treetip_left", "treetip_right"), remove = FALSE, sep = "__", fill = "right") %>%
  dplyr::mutate(`geneid` = case_when(str_detect(`treetip`, "^\\d") == TRUE ~ str_replace(`treetip`, "_.*", ""),
                                     str_detect(`treetip`, "^\\d") == FALSE ~ `treetip_left`)) %>%
  dplyr::left_join(treetip2rpsC, by = c("geneid" = "geneid")) %>%
  dplyr::left_join(bioproject2magid2domain, by = c("treetip_right" = "mag_id")) %>%
  dplyr::mutate(`treetip.new` = case_when(is.na(`domain`) == TRUE ~ paste0(`taxonid`, "_", `treetip`),
                                          `domain` != "NA" ~ paste0(`domain`, "_", `treetip`)))

tree$tip.label = tree$tip.label %>% tibble::as_tibble() %>%
  dplyr::rename(treetip = value) %>%
  dplyr::left_join(tree.tips.rename, by = c("treetip" = "treetip")) %>%
  dplyr::pull(`treetip.new`)

tips2drop = c("Archaea_QMTI01000022.1_5__GCA_003649945.1_ASM364994v1",
"Archaea_QMTF01000197.1_4__GCA_003649895.1_ASM364989v1",
"Archaea_QMVS01000135.1_1__GCA_003661225.1_ASM366122v1",
"Bacteria_QMPS01000173.1_6__GCA_003648675.1_ASM364867v1")
tree1 = di2multi(drop.tip(tree, tips2drop))
Archaeal_mrca = getMRCA(tree1, grep("Archaea", tree1$tip.label, value = T))
Bacterial_mrca = getMRCA(tree1, grep("Bacteria", tree1$tip.label, value = T))
tree2 = root(tree1, grep("Archaea", tree1$tip.label, value = T))
write.tree(tree2, file = file.path(base, "raxml_rpsC.rooted.tree"))

tree2.isolates = keep.tip(tree2, grep("^\\d", tree2$tip.label, value = T, perl = T))
write.tree(tree2.isolates, file = file.path(base, "raxml_rpsC.isolates.tree"))


# library(ggtree)
# library(phangorn)
# library(phytools)
# tree2 <- groupClade(tree1, .node = c(Archaeal_mrca, Bacterial_mrca))
#
# p = ggtree(tree1)
# p2 <- p %>% collapse(node=Archaeal_mrca) +
#   geom_point2(aes(subset=(node==Archaeal_mrca)), shape=21, size=1, fill='green')
# p2 <- collapse(p2, node=Bacterial_mrca) +
#   geom_point2(aes(subset=(node==Bacterial_mrca)), shape=23, size=5, fill='red')
# print(p2)
#
# Archaeal_mrca_desc = getDescendants(tree1, Archaeal_mrca)
# Bacterial_mrca_desc = getDescendants(tree, Bacterial_mrca)
#
# Archaeal_mrca = findMRCA(tree1, grep("Archaea", tree1$tip.label, value = T))
# Bacterial_mrca = findMRCA(tree1, grep("Bacteria", tree1$tip.label, value = T))
#
#
# library(diversitree)
# Archaeal_mrca = mrca.tipset(tree, grep("Archaea", tree$tip.label, value = T))

