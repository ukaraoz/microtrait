library(easypackages)
libraries("dplyr", "readr", "pheatmap", "corrplot", "janitor", "vegan", "factoextra", "cluster", "kmed")

calc_mixeddist = function(feature_matrix, method = "wishart", binarytype = "logical", byrow = 1) {
  # Calculate distances for mixed variable data such as
  # Gower, Podani, Wishart, Huang, Harikumar-PV, and Ahmad-Dey
  #methods = c("gower", "wishart", "podani", "huang", "harikumar", "ahmad")
  #d_gower = distmix(feature_matrix %>% dplyr::select(-1) %>% as.data.frame(),
  #                    method = "gower", idnum = idnum, idbin = idbin)
  matrix = feature_matrix %>% dplyr::select(-1) %>% as.data.frame()
  rownames(matrix) = feature_matrix %>% dplyr::pull(1)
  idnum = feature_matrix %>% dplyr::select(-1) %>% as.data.frame() %>%
    select(which(sapply(.,is.double))) %>%
    colnames() %>% match(colnames(feature_matrix)) -1
  if(binarytype == "logical") {
    idbin = feature_matrix %>% dplyr::select(-1) %>% as.data.frame() %>%
      select(which(sapply(.,is.logical))) %>%
      colnames() %>% match(colnames(feature_matrix)) -1
    distance = kmed::distmix(matrix,
                        method = method, idnum = idnum, idbin = idbin)
  }
  if(binarytype == "factor") {
    idcat = feature_matrix %>% dplyr::select(-1) %>% as.data.frame() %>%
      select(which(sapply(.,is.factor))) %>%
      colnames() %>% match(colnames(feature_matrix)) -1
    distance = kmed::distmix(matrix,
                             method = method, idnum = idnum, idcat = idcat)
  }
  distance = as.dist(distance)
  distance
}

plot_corrplot = function() {
  col3 = colorRampPalette(c('red', 'white', 'blue'))
  col5 = colorRampPalette(brewer.pal(n = 11, name = "RdBu"))

  pdf(height = 40, width = 40, file.path(base, paste0(dataset, ".alltraits.cov.ward-ordered.nolegend.pdf")))
  corrplot(cor_matrix,
           order = 'hclust', hclust.method = "ward",
           # c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid")
           #order = 'AOE', # angular order of the eigenvectors-Friendly 2002
           method = 'color',
           col = col3(16),
           col.lim = c(range(-0.6, range(cor_matrix)[2])),
           #addrect = 10, rect.lwd = 2,
           type = 'upper', diag = FALSE,
           cl.pos = "b", cl.cex = 3,
           tl.pos = "td", tl.cex = 1, tl.offset = 3, tl.col = "black"
           #type = 'full', tl.pos='n'
          )
  dev.off()
}

base = "/Users/ukaraoz/Work/microtrait/code/inst/extdata"
dataset = "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID"
genomeset_results_wmetadata = readRDS(file.path(base, paste0(dataset, ".microtraitresults.wmetadata.rds")))
genomeset_results_wmetadata_norm = trait.normalize(genomeset_results_wmetadata) %>% convert_traitdatatype(binarytype = "special")
genomeset_results_wmetadata_norm_fordist = trait.normalize(genomeset_results_wmetadata) %>% convert_traitdatatype(binarytype = "logical")
hmm_matrix = genomeset_results_wmetadata[["hmm_matrix"]]
rule_matrix = genomeset_results_wmetadata[["rule_matrix"]]
hmmandrule_matrix = hmm_matrix %>% dplyr::select(`id`:`OGT`) %>% dplyr::left_join(rule_matrix, by = c("id" = "id"))

rows = 1:nrow(genomeset_results_wmetadata_norm[[3]])
#rows = 1:4000
#glimpse(genomeset_results_wmetadata_norm[["trait_matrixatgranularity3"]])
feature_matrix = genomeset_results_wmetadata_norm[["trait_matrixatgranularity3"]] %>%
  dplyr::slice(rows) %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  dplyr::select(-`Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation`) %>%
  dplyr::select(c(1:189))

feature_matrix_fordist = genomeset_results_wmetadata_norm_fordist[["trait_matrixatgranularity3"]] %>%
  dplyr::slice(rows) %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  dplyr::select(-`Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation`) %>%
  dplyr::select(c(1:189))

#distance = calc_mixeddist(feature_matrix)

# testing
ngenomes = 1000
feature_matrix_subset = feature_matrix %>% slice(1:ngenomes)
feature_matrix_subset_dist = calc_mixeddist(feature_matrix_fordist %>% slice(1:ngenomes), binarytype = "logical")

toplot = feature_matrix_subset %>% dplyr::select(-id) %>% as.data.frame()
rownames(toplot) = feature_matrix_subset %>% dplyr::pull(id)
qnt <- quantile(as.numeric(data.matrix((toplot))), c(0.01, 0.99))
#qnt <- quantile(as.numeric(data.matrix((toplot)))[which(as.numeric(data.matrix((toplot))) >=0)], c(0.01, 0.99))
brks <- seq(qnt[1], qnt[2], length = 20)
seqPal5 <- colorRampPalette(c("black", "navyblue",
                              "mediumblue", "dodgerblue3", "aquamarine4", "green4",
                              "yellowgreen", "yellow"))(length(brks) - 1)

pheatmap(toplot, cluster_rows = TRUE, cluster_cols = TRUE,
         #treeheight_row = 0, treeheight_col = 0,
         clustering_distance_rows = feature_matrix_subset_dist,
         color = seqPal5, breaks = brks,
         height = 80, width = 50,
         fontsize_row = 6, fontsize_col = 8, angle_col = 90, border_color = "black",
         filename = file.path("~/Desktop", paste0("test", ".alltraits.heatmap.pdf")))

#ngenomes = nrow(feature_matrix)
ngenomes = 1000
toplot = feature_matrix %>% dplyr::select(-id) %>% slice(1:ngenomes) %>% as.data.frame()
#saveRDS(toplot, file.path(base, paste0(dataset, ".", "toplot.rds")))
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
