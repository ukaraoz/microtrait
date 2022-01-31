library(easypackages)
libraries("dplyr", "readr", "pheatmap", "pspearman", "RColorBrewer", "corrplot", "janitor", "vegan", "factoextra", "cluster", "kmed")

base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))

# normalize by genome size and code binary traits as factor
genomeset_results_wmetadata_norm = genomeset_results_wmetadata %>% trait.normalize
hmm_matrix = genomeset_results_wmetadata[["hmm_matrix"]]
rule_matrix = genomeset_results_wmetadata[["rule_matrix"]]
hmmandrule_matrix = hmm_matrix %>% dplyr::select(`id`:`OGT`) %>% dplyr::left_join(rule_matrix, by = c("id" = "id"))

rows = 1:nrow(genomeset_results_wmetadata_norm[[3]]) # 10305
#rows = 1:4000
#glimpse(genomeset_results_wmetadata_norm[["trait_matrixatgranularity3"]])
# filter to soil genomes, i.e. exclude aquatic
feature_matrix = genomeset_results_wmetadata_norm[["trait_matrixatgranularity3"]] %>%
  dplyr::slice(rows) %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  dplyr::select(-`Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation`) %>%
  dplyr::select(c(1:189))

##################
### calculate distances based on mixed data
##################
genomeset_results_wmetadata_norm_fordist =
  genomeset_results_wmetadata %>%
  trait.normalize %>%
  convert_traitdatatype(binarytype = "logical")
feature_matrix_fordist =
  genomeset_results_wmetadata_norm_fordist[["trait_matrixatgranularity3"]] %>%
  dplyr::slice(rows) %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  dplyr::select(-`Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation`) %>%
  dplyr::select(c(1:189))
distance = feature_matrix_fordist %>% calc_mixeddist(method = "wishart", binarytype = "logical")

feature_matrix_bin = feature_matrix %>% trait.continuous2binary
pheatmapout = cluster_traitmatrix(feature_matrix_bin,
                     clustering_distance_rows = distance,
                     clustering_distance_cols = "binary",
                     filename = file.path(base, paste0(dataset, ".alltraits.heatmap.pdf")))

hclust_rows = pheatmapout[[1]]
#which.min(abs(maov_results[,2]-70))
hcut = cutree(hclust_rows, 194)
genome2guild = data.frame(genome = names(hcut), guild = paste0("guild_", hcut))
rownames(genome2guild) = names(v)

temp = table(genome2guild[, "guild"])
guild2size = data.frame(guild = names(temp), size = as.numeric(temp)) %>% as_tibble()
genome2guild = genome2guild %>% as_tibble() %>%
  dplyr::left_join(guild2size, by = c("guild" = "guild"))
pdf(width = 20, height = 8, file = file.path(base, paste0(dataset, ".alltraits.guildsizedist.pdf")))
barplot(table(guild2size[,2]), xlab = "guild size (number of genomes)", ylab = "number of guilds", cex.lab = 1.2, cex.axis = 1.2)
dev.off()

select_guilds = guild2size %>% dplyr::filter(size>=100) %>% dplyr::pull(guild)
for(s in 1:length(select_guilds)) {
  genome2guild %>% dplyr::filter(guild == select_guilds[s])
}

feature_matrix_bin_select_guilds = feature_matrix_bin %>%
  dplyr::left_join(genome2guild, by = c("id" = "genome")) %>%
  dplyr::select(c("id", "guild", "size", everything())) %>%
  dplyr::filter(size >= 100) %>%
  dplyr::select(-c(id, size)) %>%
  dplyr::group_by(guild) %>%
  dplyr::summarise(across(everything(), mean)) %>% as.data.frame()
feature_matrix_bin_select_guilds_toplot = feature_matrix_bin_select_guilds %>% dplyr::select(-guild) %>% as.data.frame()

trait_order = pheatmapout[[2]]$labels[pheatmapout[[2]]$order]
pheatmap(feature_matrix_bin_select_guilds_toplot[, trait_order],
         color=colorRampPalette(brewer.pal(n = 7, name = "Reds"))(20),
         width = 44, height = 16,
         cellwidth = 16, cellheight = 16,
         fontsize_row = 6, fontsize_col = 14, angle_col = 90, border_color = "black",
         cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_row = c(1,2,3,4,5,6,7,8),
         show_rownames = F, show_colnames = T,
         filename = file.path(base, paste0(dataset, ".selectguild.profiles.pdf")))


toplot = feature_matrix_bin %>% dplyr::select(-id) %>% as.data.frame()
pheatmap(toplot,
         kmeans_k = 194,
         #cutree_rows = cutree_rows, cutree_cols = cutree_cols,
         clustering_distance_rows = distance,
         clustering_distance_cols = "binary",
         color = c("white", "red"), breaks = c(0,0.9,1),
         show_rownames = F, show_colnames = T,
         treeheight_row = 150, treeheight_col = 200,
         #cellwidth = 1, cellheight = 1, labels_row = NA, border_color = NA, treeheight_row = 50,
         height = 80, width = 30,
         fontsize_row = 6, fontsize_col = 8, angle_col = 90, border_color = "black",
         cluster_cols = T,
         #annotation_row = annotation_row, annotation_colors = annotation_colors,
         silent = T, legend = F,
         filename = sub(".pdf", ".kmeans.pdf", filename))
saveRDS(pheatmapout, file = "/Users/ukaraoz/Desktop/adonis_cori/pheatmapout.rds")
# compute variance across difference cuts
hclust_rows = pheatmapout[[1]]

#adonis_cori = list()
#adonis_cori[["distance"]] = distance
#adonis_cori[["pheatmapout"]] = pheatmapout
#saveRDS(adonis_cori, file = "/Users/ukaraoz/Desktop/adonis_cori/adonis_cori.rds")
# /global/u2/u/ukaraoz/bin/anaconda3/envs/R4/bin/R
library(vegan)
library(pheatmap)
library(parallel)
adonis_cori = readRDS("/global/homes/u/ukaraoz/cscratch/adonis_microtrait/adonis_cori.rds")

distance = adonis_cori[["distance"]]
pheatmapout = adonis_cori[["pheatmapout"]]
hclust_rows = pheatmapout[[1]]
nclusters = seq(2, 3556, 2)
adonis_results_all =
  parallel::mclapply(2:1778,
                     function(i) {
                       v = cutree(hclust_rows, nclusters[i])
                       genome2cluster = data.frame(cluster = factor(v))
                       rownames(genome2cluster) = names(v)
                       adonis_results = adonis2(distance ~ cluster, data = genome2cluster, perm = 1)
                       adonis_results
                     },
                    mc.cores = 25)
adonis_results1 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results1", ".rds"))
adonis_results2 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results2", ".rds"))
adonis_results3 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results3", ".rds"))
adonis_results4 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results4", ".rds"))
adonis_results5 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results5", ".rds"))
adonis_results6 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results6", ".rds"))
adonis_results7 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results7", ".rds"))
adonis_results8 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results8", ".rds"))
adonis_results9 = readRDS(paste0("/Users/ukaraoz/Work/microtrait/code/inst/extdata/adonis/adonis_results9", ".rds"))
adonis_resultslist = c(adonis_results1, adonis_results2,
                   adonis_results3, adonis_results4,
                   adonis_results5, adonis_results6,
                   adonis_results7, adonis_results8, adonis_results9)
saveRDS(adonis_resultslist, "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait-out/soilgenomes_maov_results_list.rds")
adonis_results = matrix(nrow = 0, ncol = 4)
nclusters = seq(2, 3556, 2)
for(i in 1:length(adonis_resultslist)) {
  cat(i, "\n")
  adonis_results = rbind(adonis_results,
                         c(nclusters[i], adonis_resultslist[[i]]$R2))
}

adonis_results = data.frame(`number of guilds` = adonis_results[,1],
                            `percent variance` = round(adonis_results[,2]*100,2),
                            check.names = F)
library(ggplot2)
p = ggplot(adonis_results,aes(x=`number of guilds`, y=`percent variance`)) +
  geom_line() + geom_point(size=0.8) +
  #xlim(1, 500) +
  geom_hline(yintercept = 60, color = "red", size = 0.3) +
  geom_vline(xintercept = 100, color = "red", size = 0.3) +
  geom_hline(yintercept = 70, color = "red", size = 0.3) +
  geom_vline(xintercept = 194, color = "red", size = 0.3) +
  geom_hline(yintercept = 80.04, color = "red", size = 0.3) +
  geom_vline(xintercept = 400, color = "red", size = 0.3) +
  scale_x_continuous(breaks = c(seq(0, 3556, by = 200), 3556)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))
ggsave(p, height = 6, width = 12,
       filename = file.path(base, paste0(dataset, ".adonis.percentvariance.pdf")))

{# define dendrogram object to play with:
hc <- hclust(alltraits.dist, "average")
#dend <- as.dendrogram(hc)
#results = matrix(nrow = 0, ncol =3)
cl <- parallel::makeForkCluster(30)
doParallel::registerDoParallel(cl)
results = foreach(i = 2:attr(alltraits.dist, "Size"), .combine = 'rbind') %dopar% {
  cat(i, "\n")
  o.relgr = cutree(hc, k=i)
  o.adonis = adonis(alltraits.dist ~ as.factor(o.relgr))
  c(i, o.adonis$aov.tab$SumsOfSqs[c(1,3)])
}
saveRDS(results, file = "results.jaccarddist.adonis.rds")
parallel::stopCluster(cl)

for(i in 2:attr(alltraits.dist, "Size")){
  cat(k, "\n")
  o.relgr = cutree(hc, k=i)
  o.adonis = adonis(alltraits.dist ~ as.factor(o.relgr))
  #o.adonis$aov.tab$SumsOfSqs[c(1,3)]
  results = rbind(results,
                  c(k, o.adonis$aov.tab$SumsOfSqs[c(1,3)]))
}

for(i in 1:9) {
  temp = readRDS(paste0("adonis_results", i, ".rds"))
}
}
