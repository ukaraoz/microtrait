library(dplyr)
library(ggplot2)
library(RColorBrewer)

base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata/traitconservation"
models = readRDS(file.path(base, "nlme.models.optimizer-nlm.species.rds"))
#models1 = readRDS(file.path(base, "nlme.models.rds"))
#nlme.models.optimizer-nlm.rds

traits_wnospace = c("mingentime", "OGT", traits_listbygranularity[[3]] %>%
             #dplyr::filter(`microtrait_trait-type` %in% c("count","count_by_substrate")) %>%
             dplyr::select(`microtrait_trait-name`) %>% pull() %>% as.character()) %>% tibble::as_tibble() %>%
  dplyr::rename(`microtrait_trait-name` = `value`) %>%
  dplyr::mutate(`microtrait_trait-name-nospace` = gsub(" |:|-|,|/|\\(|\\)", "_", `microtrait_trait-name`))
write.table(traits_wnospace %>% as.data.frame(),
            file = file.path(base, "traits_wnospace.xls"),
            quote = F, row.names = F, col.names = T, sep = "\t")

# remove traits with errors
indices.select = c()
for(i in 1:length(models)) {
  temp = models[[i]] %>% dplyr::filter(`variable` == "error") %>% dplyr::select(`variable`)
  if(nrow(temp) == 0) {
    indices.select = c(indices.select, i)
  }
}

# a temporary fix
for(i in 1:length(indices.select)) {
  model1_within=100-models[[indices.select[i]]] %>% dplyr::filter(`model`==model1) %>% dplyr::select(variance) %>% sum()
  toadd1 = data.frame(variable = "Within",
                      variance = model1_within,
                      trait = models[[indices.select[i]]] %>% dplyr::pull(`trait`) %>% unique(),
                      model = model1)

  model2_within=100-models[[indices.select[i]]] %>% dplyr::filter(`model`==model2) %>% dplyr::select(variance) %>% sum()
  toadd2 = data.frame(variable = "Within",
                      variance = model2_within,
                      trait = models[[indices.select[i]]] %>% dplyr::pull(`trait`) %>% unique(),
                      model = model2)
  models[[indices.select[i]]] = dplyr::bind_rows(models[[indices.select[i]]], toadd1, toadd2)
}
model1 = "~1 | NCBI_Phylum/NCBI_Class/NCBI_Order/NCBI_Family/NCBI_Genus/NCBI_Species/Ecosystem_Category"
model2 = "~1 | Ecosystem_Category/NCBI_Phylum/NCBI_Class/NCBI_Order/NCBI_Family/NCBI_Genus/NCBI_Species"

# create a dataset
toplot = do.call(bind_rows, models[indices.select])
toplot.model1 = toplot %>%
  dplyr::filter(model == model1)
traits.ordered = toplot.model1 %>%
  dplyr::filter(`variable` == "Within") %>%
  dplyr::arrange(variance) %>%
  dplyr::select(`trait`) %>%
  dplyr::left_join(traits_wnospace, by = c("trait" = "microtrait_trait-name-nospace"))

toplot.model1 = toplot.model1 %>%
  dplyr::mutate(trait = factor(toplot.model1$trait,
                               levels = traits.ordered %>% dplyr::pull(`trait`),
                               labels = traits.ordered %>% dplyr::pull(`microtrait_trait-name`),
                               ordered = T),
                variable = factor(toplot.model1$variable,
                                  levels = c("NCBI_Phylum", "NCBI_Class", "NCBI_Order", "NCBI_Family", "NCBI_Genus", "NCBI_Species", "Ecosystem_Category", "Within"),
                                  labels = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Ecosystem", "Unexplained"),
                                  ordered = T),
                variance = as.numeric(toplot.model1$variance)
                )

# Stacked + percent
# brewer.pal(9, "YlGnBu")
p = ggplot(toplot.model1, aes(fill=variable, width=.9)) +
  geom_col(aes(x=trait, y=variance)) +
  #scale_fill_manual(values=c("#081D58","#253494","#225EA8","#1D91C0","#41B6C4","red","gray")) +
  scale_fill_manual(values=c("#17406D", "#0F6FC6", "#009DD9", "#0BD0D9", "#10CF9B", "#7CCA62","red","gray")) +
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5, hjust=1))
ggsave(p, width = 26, height = 16, file = file.path(base, "nlme.model1.pdf"))

toplot.model2 = toplot %>%
  dplyr::filter(model == model2)
traits.ordered = toplot.model2 %>%
  dplyr::filter(`variable` == "Within") %>%
  dplyr::arrange(variance) %>%
  dplyr::select(`trait`) %>%
  dplyr::left_join(traits_wnospace, by = c("trait" = "microtrait_trait-name-nospace"))

toplot.model2 = toplot.model2 %>%
  dplyr::mutate(trait = factor(toplot.model2$trait,
                               levels = traits.ordered %>% dplyr::pull(`trait`),
                               labels = traits.ordered %>% dplyr::pull(`microtrait_trait-name`),
                               ordered = T),
                variable = factor(toplot.model2$variable,
                                  levels = c("NCBI_Phylum", "NCBI_Class", "NCBI_Order", "NCBI_Family", "NCBI_Genus", "NCBI_Species", "Ecosystem_Category", "Within"),
                                  labels = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Ecosystem", "Unexplained"),
                                  ordered = T),
                variance = as.numeric(toplot.model2$variance)
  )
# Stacked + percent
# brewer.pal(9, "YlGnBu")
p = ggplot(toplot.model2, aes(fill=variable, width=.9)) +
  geom_col(aes(x=trait, y=variance)) +
  #scale_fill_manual(values=c("#081D58","#253494","#225EA8","#1D91C0","#41B6C4","red","gray")) +
  scale_fill_manual(values=c("#17406D", "#0F6FC6", "#009DD9", "#0BD0D9", "#10CF9B", "#7CCA62","red","gray")) +
                             theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5, hjust=1))
ggsave(p, width = 26, height = 16, file = file.path(base, "nlme.model2.pdf"))

towrite1 = toplot.model1 %>%
  tidyr::pivot_wider(names_from = c("variable","model"),
                     values_from = "variance",
                     values_fill = 0) %>%
  dplyr::mutate(trait = levels(trait)[as.numeric(trait)])
towrite2 = toplot.model2 %>%
  tidyr::pivot_wider(names_from = c("variable","model"),
                     values_from = "variance",
                     values_fill = 0) %>%
  dplyr::mutate(trait = levels(trait)[as.numeric(trait)])
towrite = towrite1 %>% dplyr::left_join(towrite2, by = c("trait" = "trait"))
write.table(towrite %>% as.data.frame(),
            file = file.path(base, "traitconservation_bytaxonomy.xls"),
            quote = F, row.names = F, col.names = T, sep = "\t")
