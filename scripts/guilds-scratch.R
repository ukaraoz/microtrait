for(m in 1:5) {
  # feature_matrix %>%
  #   select_if(function(x) any(is.na(x))) %>%
  #   summarise_each(funs(sum(is.na(.)))) -> extra_NA
  cat(m, "\n")
  d = distmix(feature_matrix %>% dplyr::select(-1) %>% as.data.frame(),
              method = methods[m], idnum = idnum, idbin = idbin)
  assign(paste0("d_", methods[m]), d)
  # gower_dist <- cluster::daisy(feature_matrix[, -1],
  #                              metric = "gower")
}

ade4::mantel.rtest(as.dist(d_gower), as.dist(d_gower))
vegan::mantel(as.dist(d_gower), as.dist(d_wishart))

vegan::mantel(as.dist(d_gower[1:1000,1:1000]), as.dist(d_gower[1:1000,1:1000]))
#toplot = feature_matrix
feature_matrix_df = feature_matrix %>% dplyr::select(-1) %>% as.data.frame()

feature_matrix_bin = trait.continuous2binary(feature_matrix)
feature_matrix_bin_df = feature_matrix_bin %>% dplyr::select(-id) %>% as.data.frame()
range = 1:ncol(feature_matrix_bin_df)
cor_matrix = cor(feature_matrix_bin_df[, range])
cor_melted = melt_dist(cor_matrix)
pdf(height = 6, width = 12,
    file.path(base, paste0(dataset, ".", "trait_matrixatgranularity3", "_traitcorr-hist.pdf")))
h = hist(cor_melted[,"dist"], breaks = seq(-0.6, 1, 0.1), xlab = "trait-trait correlation", ylab = "count", main = "", cex.lab = 1.5,xaxt= "n")
axis(side=1, at = c(seq(-0.6, 0, 0.1), seq(0.1, 1, 0.1)), labels = c(seq(-0.6, 0, 0.1), seq(0.1, 1, 0.1)))
dev.off()
write.table(cor_melted,
            file = file.path(base, paste0(dataset, ".", "trait_matrixatgranularity3", "_traitcorr.xls")),
            row.names = F, col.names = T, sep = "\t", quote = F)

# cluster soil genomes
# Aquatic: 4698, Terrestrial: 3693, Plants: 376
feature_matrix = genomeset_results_wmetadata[["trait_matrixatgranularity3"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  dplyr::select(-`Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation`) %>%
  dplyr::select(c(1:189))
feature_matrix_bin = trait.continuous2binary(feature_matrix)
feature_matrix_bin_df = feature_matrix_bin %>% dplyr::select(-id) %>% as.data.frame()
saveRDS(feature_matrix_bin_df, file.path(base, paste0(dataset, ".", "feature_matrix_bin_soil_df.rds")))

compute.prevalence.help

ngenomes = nrow(feature_matrix_bin)
#ngenomes = 200
toplot = feature_matrix_bin %>% dplyr::select(-id) %>% slice(1:ngenomes) %>% as.data.frame()
saveRDS(toplot, file.path(base, paste0(dataset, ".", "toplot.rds")))
pheatmap.outlist = list()
# all: k=277; 70.19208
# soil: k=192; 70.01773
pheatmapout = pheatmap(toplot,
                       #kmeans_k = k, iter.max = 10, nstart = 25,
                       cutree_rows = 192, cutree_cols=15,
                       color = c("white", "red"), breaks = c(0,0.9,1),
                       show_rownames = F, show_colnames = T,
                       treeheight_row = 150, treeheight_col = 200,
                       #cellwidth = 1, cellheight = 1, labels_row = NA, border_color = NA, treeheight_row = 50,
                       height = 80, width = 30,
                       fontsize_row = 6, fontsize_col = 8, angle_col = 90, border_color = "black",
                       cluster_cols = T, clustering_distance_rows = "binary",
                       #annotation_row = annotation_row, annotation_colors = annotation_colors,
                       silent = T, legend = F,
                       filename = file.path(base, paste0(dataset, ".alltraits.heatmap.pdf"))
)

hc = hclust(dist(toplot, method = "binary"), method = "complete")
v = cutree(hc, 192)[hc$order]
table(v)

cutree_rows
#c(totss = pheatmapout$kmeans$totss,
#  tot.withinss = pheatmapout$kmeans$tot.withinss,
#  betweenss = pheatmapout$kmeans$betweenss)
#pheatmapout$kmeans$cluster



for(k in 2:ngenomes) {
  cat(k, "\n")
  temp = kmeans(toplot, k, iter.max = 10, nstart = 25)
  between_over_tot = temp$betweenss/temp$totss*100

  pheatmap.outlist[[as.character(k)]][["ss"]] = c(totss = temp$totss,
                                                  tot.withinss = temp$tot.withinss,
                                                  betweenss = temp$betweenss)
  pheatmap.outlist[[as.character(k)]][["cluster"]] = temp$kmeans$cluster
}





p = factoextra::fviz_nbclust(toplot, kmeans, method = "wss", k.max = 100)
cluster::clusGap()

between_SS / total_SS

# gap_stat <- cluster::clusGap(toplot, FUN = kmeans, iter.max = 10, nstart = 25,
#                     K.max = 200, B = 20)
# factoextra::fviz_gap_stat(gap_stat)




temp = feature_matrix_bin_df[1:200,]
d = vegan::vegdist(temp)




M1 = cor(mtcars)
corrplot(M1, method = 'color', diag = FALSE, type = 'upper', tl.cex = 9) # colorful number


detach("package:pheatmap", unload=TRUE)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  #value = "draw_colnames",
  ns = asNamespace("pheatmap")
)

pheatmap(
  mat               = mat,
  color             = inferno(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  cluster_cols      = mat_cluster_cols,
  cluster_rows      = mat_cluster_rows,
  cellwidth         = 20,
  show_colnames     = TRUE,
  show_rownames     = FALSE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Rotated Column Names"
)


alltraits.dist = vegan::vegdist(feature_matrix, method = "jaccard")

toplot1 = toplot %>%
  dplyr::mutate_at(vars(!starts_with("id")),
                   #funs(case_when(. >= 1 ~ 1,TRUE ~ 0))) %>%
                   .funs = list(binary = ~case_when(. >= 22 ~ 1,TRUE ~ 0)))

compute.prevalence.help <- function(data) {
  prevalence = data %>% tibble::as_tibble() %>%
    # dplyr::select(c("id",
    #                 traits_listbygranularity[[3]] %>%
    #                   #dplyr::filter(`microtrait_trait-type` == "binary") %>%
    #                   dplyr::pull(`microtrait_trait-name`))) %>%
    dplyr::mutate_at(vars(!starts_with("id")),
                     funs(case_when(. >= 1 ~ 1,TRUE ~ 0))) %>%
    tidyr::pivot_longer(!id, names_to = "trait", values_to = "presence") %>%
    dplyr::select(-id) %>%
    dplyr::count(trait, presence,.drop=FALSE) %>%
    dplyr::group_by(trait) %>%  #methaneoxidation_genes = c("mmoX", "mmoY", "mmoZ", "mmoC", "mmoD", "pmoA-amoA", "pmoB-amoB", "pmoC-amoC")
    dplyr::mutate(percent = n/sum(n)*100) %>%
    dplyr::filter(`presence` == 1) %>%
    # dplyr::mutate(`percent` = 100-`percent`) %>%
    dplyr::select(c("trait", "n", "percent"))
  prevalence
}


#####cori
toplot = readRDS("organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID.feature_matrix_bin_soil_df.rds")
compute_variance_results = parallel::mclapply(2:300,
                                              function(i) {
                                                returnList = compute_variance(toplot, i)
                                                returnList
                                              },
                                              mc.cores = 56)
#compute_variance_results[[50]][[1]]$ss
compute_variance = function(data, k) {
  temp = kmeans(toplot, k, iter.max = 10, nstart = 25)
  result = list()
  result[[as.character(k)]][["ss"]] = c(totss = temp$totss,
                                        tot.withinss = temp$tot.withinss,
                                        betweenss = temp$betweenss,
                                        betweenovertot_percent = temp$betweenss/temp$totss*100
  )
  result[[as.character(k)]][["cluster"]] = temp$cluster
  result
}
#####cori
