library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(foreach)
library(doParallel)
library(pals)
library(RColorBrewer)
library(pheatmap)

base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))
objects = c("trait_matrixatgranularity3", "trait_matrixatgranularity2", "trait_matrixatgranularity1", "rule_matrix")
i = 3
# Gianna: amino acids, organic acids, nucleotides, sugars, auxins, fatty acids
select_traits_granularity1 = c("Resource Acquisition:Substrate uptake:free amino acids transport",
"Resource Acquisition:Substrate uptake:carboxylate transport",
"Resource Acquisition:Substrate uptake:nucleic acid component transport",
"Resource Acquisition:Substrate uptake:carbohydrate transport",
"Resource Acquisition:Substrate uptake:lipid transport",
"Resource Acquisition:Substrate uptake:aromatic acid transport",
"Resource Acquisition:Substrate uptake:peptide transport")
select_traits_granularity3 = c("Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:cellulose breakdown",
"Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:chitin breakdown",
"Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:heteromannan breakdown",
"Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:mixed linkage glucan breakdown",
"Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:xylan and heteroxylan breakdown",
"Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:xyloglucan breakdown")
MGT_OGT = c("mingentime", "OGT")
select_traits = c(select_traits_granularity1, select_traits_granularity3)
select_metadata = c("Ecosystem", "Ecosystem_Category", "Ecosystem_Type", "Ecosystem_Subtype", "Specific Ecosystem",
                    "NCBI_Superkingdom", "NCBI_Kingdom", "NCBI_Phylum", "NCBI_Class",
                    "NCBI_Order", "NCBI_Family", "NCBI_Genus", "NCBI_Species", "Estimated Size")

feature_matrix_matrixatgranularity1 = genomeset_results_wmetadata[["trait_matrixatgranularity1"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Aquatic", "Terrestrial", "Plants")) %>%
  dplyr::select(c("id", select_traits_granularity1))
feature_matrix_matrixatgranularity3 = genomeset_results_wmetadata[["trait_matrixatgranularity3"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Aquatic", "Terrestrial", "Plants")) %>%
  dplyr::select(c("id", select_traits_granularity3, MGT_OGT, select_metadata))
feature_matrix = feature_matrix_matrixatgranularity1 %>%
  dplyr::inner_join(feature_matrix_matrixatgranularity3, by = c("id" = "id"))
habitats_temp = feature_matrix %>% dplyr::select("Ecosystem_Category", "Ecosystem_Type") %>%
  dplyr::mutate(Ecosystem_Type_long = paste0(`Ecosystem_Category`, ":", `Ecosystem_Type`))
toplot = data.frame(Ecosystem = names(table(habitats_temp[, "Ecosystem_Type_long"])),
                    Genomes = as.numeric(table(habitats_temp[, "Ecosystem_Type_long"])))

p<-ggplot(data=habitats_temp, aes(Ecosystem_Type_long)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#trait_cols = 2

#select_traits = c(grep("transporter", colnames(taxonids2traitwmetadata), value = T),
#                  "xylanandheteroxylan", "cellulose", "xyloglucan", "chitin", "heteromannan")
# normalized by genome size
feature_matrix_norm = feature_matrix %>% as.data.frame()
for(r in 1:nrow(feature_matrix_norm)) {
  cat(r, "\n")
  feature_matrix_norm[r, select_traits] =
    (feature_matrix_norm[r, select_traits]/feature_matrix_norm[r,"Estimated Size"])*1E8
}
feature_matrix_norm = feature_matrix_norm %>% tibble::as_tibble()
feature_matrix_norm = feature_matrix_norm %>%
  dplyr::rowwise() %>%
  dplyr::mutate(sum = sum(c_across(`Resource Acquisition:Substrate uptake:free amino acids transport`:`Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:xyloglucan breakdown`))) %>%
  dplyr::filter(sum != 0) %>%
  dplyr::ungroup()

# z-score transform
feature_matrix_norm_z = feature_matrix_norm %>%
  dplyr::mutate_at(select_traits, scale)

# show genomic investment distributions
temp = feature_matrix_norm %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Marine", "Freshwater", "Soil", "Rhizoplane", "Rhizosphere", "Sediment")) %>%
  dplyr::select(c(1:14, 17)) %>%
  tidyr::pivot_longer(-c(`id`, `Ecosystem_Type`), names_to = "trait", values_to = "count per bp.") %>%
  dplyr::mutate(trait = factor(trait, levels = select_traits, labels = sub(".*:", "", select_traits), ordered = T))
p1 = ggplot(temp, aes(x=Ecosystem_Type, y=`count per bp.`, group = Ecosystem_Type)) +
  geom_boxplot() + facet_grid(. ~ trait) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.x = element_text(hjust=0, size = 10)) +
  ylab("count per bp. x 1E8 (a.u.)")
p2 = ggplot(temp, aes(x=trait, y=`count per bp.`, group = trait)) +
  geom_boxplot() + facet_grid(. ~ Ecosystem_Type) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.x = element_text(hjust=0, size = 10)) +
  ylab("count per bp. x 1E8 (a.u.)")
ggsave(p1, width = 20, height = 8,
       file = "/Users/ukaraoz/Work/microtrait/code/inst/extdata/tradeoffs/resourceacquisition1.pdf")
ggsave(p2, width = 12, height = 8,
       file = "/Users/ukaraoz/Work/microtrait/code/inst/extdata/tradeoffs/resourceacquisition2.pdf")
# categorize/bin
feature_matrix_norm_z_binned <- feature_matrix_norm_z %>%
  mutate(`Resource Acquisition:Substrate uptake:free amino acids transport binned`= cut(x = `Resource Acquisition:Substrate uptake:free amino acids transport`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate uptake:carboxylate transport binned`= cut(x = `Resource Acquisition:Substrate uptake:carboxylate transport`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate uptake:nucleic acid component transport binned`= cut(x = `Resource Acquisition:Substrate uptake:nucleic acid component transport`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate uptake:carbohydrate transport binned`= cut(x = `Resource Acquisition:Substrate uptake:carbohydrate transport`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate uptake:lipid transport binned`= cut(x = `Resource Acquisition:Substrate uptake:lipid transport`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate uptake:aromatic acid transport binned`= cut(x = `Resource Acquisition:Substrate uptake:aromatic acid transport`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate uptake:peptide transport binned`= cut(x = `Resource Acquisition:Substrate uptake:peptide transport`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:cellulose breakdown binned`= cut(x = `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:cellulose breakdown`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:chitin breakdown binned`= cut(x = `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:chitin breakdown`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:heteromannan breakdown binned`= cut(x = `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:heteromannan breakdown`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:mixed linkage glucan breakdown binned`= cut(x = `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:mixed linkage glucan breakdown`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:xylan and heteroxylan breakdown binned`= cut(x = `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:xylan and heteroxylan breakdown`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")),
         `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:xyloglucan breakdown binned`= cut(x = `Resource Acquisition:Substrate degradation:complex carbohydrate depolymerization:xyloglucan breakdown`,
                 breaks=c(-Inf, -1.5, 1.5, Inf), labels=c("-1","0","1")))
feature_matrix_norm_z_binned = mutate_if(feature_matrix_norm_z_binned, is.factor, ~ as.integer(as.character(.x)))
feature_matrix_norm_z_binned_soil = feature_matrix_norm_z_binned %>%
  dplyr::filter(Ecosystem_Type == "Soil")


# Display row and color annotations
pheatmap(test, annotation_col = annotation_col)
pheatmap(test, annotation_col = annotation_col, annotation_legend = FALSE)
pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row)

rownames(feature_matrix_norm_z_binned_soil) = feature_matrix_norm_z_binned_soil[, "id"]
pheatmap(mat = feature_matrix_norm_z_binned_soil[, paste0(select_traits, " binned")],
         #color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))
         color = c("blue", "white", "red"),
         show_rownames = F, show_colnames = F,
         cluster_cols = F,
         #cellwidth = 1, cellheight = 1,
         width = 8, height = 80, fontsize = 8, angle_col = 90, border = "black",
         filename = "/Users/ukaraoz/Work/microtrait/code/inst/extdata/tradeoffs/transportervsextracellular_enrichment_all.pdf")

pheatmap_df = feature_matrix_norm_z_binned_soil[which(apply(feature_matrix_norm_z_binned_soil[,paste0(select_traits, " binned")],1,sum)!=0), c("id", paste0(select_traits, " binned"), "mingentime", "OGT", "Estimated Size")] %>% as.data.frame(check.names = F)
rownames(pheatmap_df) = pheatmap_df[, "id"]
# row annotations, MGT and OGT
annotation_row = data.frame(feature_matrix_norm_z_binned_soil[, c("mingentime", "OGT", "Estimated Size")])
rownames(annotation_row) = feature_matrix_norm_z_binned_soil %>% pull(`id`)
colnames(annotation_row)[3] = "Genome Size"
pheatmap(mat = pheatmap_df[, paste0(select_traits, " binned")],
         #color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))
         color = c("blue", "white", "red"),
         show_rownames = F, show_colnames = F,
         cluster_cols = F,
         #cutree_rows = 10,
         #kmeans_k = 10,
         annotation_row = annotation_row[, c("mingentime", "OGT", "Genome Size")],
         #cellwidth = 12, cellheight = 12,
         width = 8, height = 12, fontsize = 8, angle_col = 90, border = "black",
         filename = "/Users/ukaraoz/Work/microtrait/code/inst/extdata/tradeoffs/transportervsextracellular_enrichment.pdf")

# manually copied cluster members, read
exudateusers = as.character(read.table("/Users/ukaraoz/Work/microtrait/code/inst/extdata/tradeoffs/exudateusers.ids.txt",sep = "\t")[,1])
exudateusers1 = as.character(read.table("/Users/ukaraoz/Work/microtrait/code/inst/extdata/tradeoffs/exudateusers.ids1.txt",sep = "\t")[,1])
polymerusers = as.character(read.table("/Users/ukaraoz/Work/microtrait/code/inst/extdata/tradeoffs/polymerusers.ids.txt",sep = "\t")[,1])
pheatmap_df = data.frame(pheatmap_df, guild = rep("", nrow(pheatmap_df)))
pheatmap_df[exudateusers, "guild"] = "exudate users"
pheatmap_df[polymerusers, "guild"] = "polymer users"
toplot = pheatmap_df[which(pheatmap_df$guild != ""), ]
ggplot(toplot, aes(x=guild, y=mingentime)) +
  geom_boxplot()
ggplot(toplot, aes(x=guild, y=`OGT`)) +
  geom_boxplot()
ggplot(toplot, aes(x=guild, y=`Estimated.Size`)) +
  geom_boxplot()

pheatmap(mat = feature_matrix_norm_z_binned_soil[which(apply(feature_matrix_norm_z_binned_soil[,paste0(select_traits, " binned")],1,sum)!=0), paste0(select_traits, " binned")],
         #color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))
         color = c("blue", "white", "red"),
         show_rownames = F, show_colnames = T,
         cluster_cols = F,
         #cellwidth = 12, cellheight = 12,
         width = 8, height = 12, fontsize = 8, angle_col = 90, border = "black",
         filename = "/Users/ukaraoz/Work/microtrait/code/inst/extdata/tradeoffs/transportervsextracellular_enrichment.wcolnames.pdf")

