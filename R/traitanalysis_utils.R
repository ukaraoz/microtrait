
cor.mixed = function(df, method.numeric = "pearson", method.logical = "cramerV") {
  # for testing
  # numeric_pos = which(sapply(df, class) == "numeric")[1:2]
  # logical_pos = which(sapply(df, class) == "logical")[1:2]
  # pos_1 = numeric_pos[1];pos_2 = numeric_pos[2];
  # r <- stats::cor(df[[pos_1]], df[[pos_2]],
  #                 method = method.numeric)
  # pos_1 = numeric_pos[1];pos_2 = logical_pos[2];
  # r <- sqrt(summary(stats::lm(df[[pos_1]] ~ as.factor(df[[pos_2]])))[["r.squared"]])
  #
  # pos_1 = logical_pos[1];pos_2 = logical_pos[2];
  # r <- lsr::cramersV(df[[pos_1]], df[[pos_2]], simulate.p.value = TRUE)

  stopifnot(inherits(df, "data.frame"))
  stopifnot(sapply(df, class) %in% c("integer",
                                     "numeric",
                                     "logical"))

  cor_fun <- function(pos_1, pos_2){
    # both are numeric
    if(class(df[[pos_1]]) %in% c("integer", "numeric") &&
       class(df[[pos_2]]) %in% c("integer", "numeric")){
      r <- stats::cor(df[[pos_1]], df[[pos_2]],
                      method = method.numeric)
    }

    # one is numeric and other is a factor/character
    if(class(df[[pos_1]]) %in% c("integer", "numeric") &&
       class(df[[pos_2]]) %in% c("logical")){
      if(nlevels(as.factor(df[[pos_2]])) == 1) {
        r = NaN
      } else {
        r <- sqrt(
          summary(
            stats::lm(df[[pos_1]] ~ as.factor(df[[pos_2]])))[["r.squared"]])
      }
    }

    if(class(df[[pos_1]]) %in% c("logical") &&
       class(df[[pos_2]]) %in% c("integer", "numeric")){
      if(nlevels(as.factor(df[[pos_1]])) == 1) {
        r = NaN
      } else {
        r <- sqrt(
          summary(
            stats::lm(df[[pos_2]] ~ as.factor(df[[pos_1]])))[["r.squared"]])
      }
    }

    # both are logical
    if(class(df[[pos_1]]) %in% c("logical") &&
       class(df[[pos_2]]) %in% c("logical")){
      if(nlevels(as.factor(df[[pos_1]])) == 1 |
         nlevels(as.factor(df[[pos_2]])) == 1) {
        r = NaN
      } else {
        if(method.logical == "cramerV") {
          r <- cramersV(df[[pos_1]], df[[pos_2]], simulate.p.value = TRUE)
        }
        if(method.logical == "phicoef") {
          r <- phicoef(df[[pos_1]], df[[pos_2]])
        }
        if(method.logical == "jaccard") {
          r <- jaccard(df[[pos_1]], df[[pos_2]])
        }
      }
    }

    return(r)
  }

  cor_fun <- Vectorize(cor_fun)

  # now compute corr matrix
  corrmat <- outer(1:ncol(df),
                   1:ncol(df),
                   function(x, y) cor_fun(x, y)
  )

  rownames(corrmat) <- colnames(df)
  colnames(corrmat) <- colnames(df)

  return(corrmat)
}

phicoef = function(x, y = NULL) {
  if (is.null(y)) {
    if (!is.integer(x) || length(x) != 4L)
      stop("when 'y' is not supplied, 'x' must be ", "a 2x2 integer matrix or an integer vector of length 4")
    a <- x[1L]
    c <- x[2L]
    b <- x[3L]
    d <- x[4L]
  }
  else {
    if (!is.logical(x) || !is.logical(y) || length(x) !=
        length(y))
      stop("when 'y' is supplied, 'x' and 'y' must be ",
           "2 logical vectors of the same length")
    a <- sum(x & y)
    b <- sum(x & !y)
    c <- sum(!x & y)
    d <- sum(!x & !y)
  }
  a <- as.double(a)
  b <- as.double(b)
  c <- as.double(c)
  d <- as.double(d)
  div <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
  (a * d - b * c)/div
}

jaccard = function(x,y) {
  d=sum(x==TRUE & y==TRUE, na.rm=TRUE)
  b=sum(x==TRUE & y==FALSE, na.rm=TRUE)
  c=sum(x==FALSE & y==TRUE, na.rm=TRUE)
  res=d/(b+c+d)
  return(res)
}

cramersV = function(...) {
  test <- stats::chisq.test(...)
  chi2 <- test$statistic
  N <- sum(test$observed)

  if (test$method =="Chi-squared test for given probabilities"){
    # for GOF test, calculate max chi-square value
    ind <- which.min(test$expected)
    max.dev <- test$expected
    max.dev[ind] <- N-max.dev[ind]
    max.chi2 <- sum( max.dev ^2 / test$expected )
    V <- sqrt( chi2 / max.chi2 )
  }
  else {
    # for test of association, use analytic expression
    k <- min(dim(test$observed))
    V <- sqrt( chi2 / (N*(k-1)) )
  }
  names(V) <- NULL
  return(V)
}

#' Calculate distances between genomes.
#'
#' @param trait_matrix_bin
#'
#' @import kmed
#' @export
calc_mixeddist = function(trait_matrix, idcol = "id", col2ignore = c("genome_length", "mingentime"), method = "wishart", binarytype = "logical", byrow = 1, verbose = TRUE) {
  # Calculate distances for mixed variable data such as
  # Gower, Podani, Wishart, Huang, Harikumar-PV, and Ahmad-Dey
  #methods = c("gower", "wishart", "podani", "huang", "harikumar", "ahmad")
  #d_gower = distmix(trait_matrix %>% dplyr::select(-1) %>% as.data.frame(),
  #                    method = "gower", idnum = idnum, idbin = idbin)
  if(verbose) {message("Calculating trait distance matrix for ", nrow(trait_matrix), " genomes.")}
  matrix = trait_matrix %>% dplyr::select(-idcol) %>% dplyr::select(-any_of(col2ignore)) %>% as.data.frame()
  rownames(matrix) = trait_matrix %>% dplyr::pull(1)
  idnum = trait_matrix %>% dplyr::select(-idcol) %>% dplyr::select(-any_of(col2ignore)) %>% as.data.frame() %>%
    select(which(sapply(.,is.double))) %>%
    colnames() %>% match(colnames(trait_matrix)) -1
  if(binarytype == "logical") {
    idbin = trait_matrix %>% dplyr::select(-idcol) %>% dplyr::select(-any_of(col2ignore)) %>% as.data.frame() %>%
      select(which(sapply(.,is.logical))) %>%
      colnames() %>% match(colnames(trait_matrix)) -1
    # distmix returns type matrix
    distance = kmed::distmix(matrix,
                             method = method, idnum = idnum, idbin = idbin)
  }
  if(binarytype == "factor") {
    idcat = trait_matrix %>% dplyr::select(-idcol) %>% dplyr::select(-any_of(col2ignore)) %>% as.data.frame() %>%
      select(which(sapply(.,is.factor))) %>%
      colnames() %>% match(colnames(trait_matrix)) -1
    distance = kmed::distmix(matrix,
                             method = method, idnum = idnum, idcat = idcat)
  }
  distance = as.dist(distance)
  distance
}

#' Analyze trait2trait correlations.
#'
#' @param trait_matrix_bin
#'
#' @import RColorBrewer corrplot ggplot2
#' @export
trait2trait_corr = function(trait_matrix_bin, idcol = "id", annot_cols = c("mingentime", "ogt"), verbose = T, outdir, dataset) {
  ##############
  ### trait-to-trait correlation
  ##############
  trait_matrix_bin_df = trait_matrix_bin %>% dplyr::select(-all_of(idcol)) %>% as.data.frame()

  # trait_matrix_bin_df = trait_matrix_bin %>% dplyr::select(-all_of(idcol)) %>% dplyr::select(-all_of(annot_cols)) %>% as.data.frame()

  if(verbose) {message}
  # remove all absent/all present traits
  ngenomes = nrow(trait_matrix_bin_df)
  ntraits = ncol(trait_matrix_bin_df)
  if(verbose) {message("Read a trait matrix with ", ngenomes, " genomes and ", ntraits, " traits dimensions.")}
  traitprevalence = apply(trait_matrix_bin_df, 2, sum)
  undetectedtraits = which(traitprevalence == 0)
  alldetectedtraits = which(traitprevalence == ngenomes)
  range = setdiff(1:ntraits, c(undetectedtraits, alldetectedtraits))
  if(verbose) {message("Removing ", length(c(undetectedtraits, alldetectedtraits)) ," undetected and uninformative traits.")}
  if(verbose) {message("Analyzing a trait matrix with ", ngenomes, " genomes and ", length(range), " traits dimensions.")}

  cor_matrix = cor(trait_matrix_bin_df[, range], method = "spearman")
  cor_melted = melt_dist(cor_matrix)
  pdf(height = 6, width = 14,
      file.path(outdir, paste0(dataset, "_traitcorr-hist.pdf")))
  h = hist(cor_melted[,"dist"], breaks = seq(-0.8, 1, 0.1), xlab = "", ylab = "", main = "", cex.axis = 1.5, las = 1, xaxt= "n")
  axis(side=1, at = c(seq(-0.6, 0, 0.1), seq(0.1, 1, 0.1)), labels = c(seq(-0.6, 0, 0.1), seq(0.1, 1, 0.1)), cex.axis = 1.5)
  dev.off()
  #if(verbose) {message("Generated trait correlation histogram: ", file.path(outdir, paste0(dataset, "_traitcorr-hist.pdf")))}

  write.table(cor_melted,
              file = file.path(outdir, paste0(dataset, "_traitcorr.xls")),
              row.names = F, col.names = T, sep = "\t", quote = F)
  if(verbose) {message("Wrote trait correlations: ", file.path(outdir, paste0(dataset, "_traitcorr.xls")))}


  col3 = colorRampPalette(c('red', 'white', 'blue'))
  col5 = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))

  pdf(height = 40, width = 40,
      file.path(outdir, paste0(dataset, "_traitcorr-plot.pdf")))
  corrplot::corrplot(cor_matrix,
           order = 'hclust', hclust.method = "ward.D",
           # c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid")
           #order = 'AOE', # angular order of the eigenvectors-Friendly 2002
           method = 'color',
           col = col3(16),
           col.lim = c(range(-0.8, range(cor_matrix)[2])),
           #addrect = 10, rect.lwd = 2,
           type = 'upper', diag = FALSE,
           cl.pos = "b", cl.cex = 3,
           tl.pos = "td", tl.cex = 1, tl.offset = 3, tl.col = "black",
           mar = c(2, 2, 2, 2)
           #type = 'full', tl.pos='n'
  )
  dev.off()
  if(verbose) {message("Generated trait correlation heatmap: ",
                       file.path(outdir, paste0(dataset, "_traitcorr-plot.pdf")))}
  #cor_melted1 = cor_melted %>% as_tibble() %>%
  #  tidyr::separate(variable1, c("strategy1"), sep = ":", remove = F) %>%
  #  tidyr::separate(variable2, c("strategy2"), sep = ":", remove = F) %>%
  #  dplyr::mutate(`group` =
  #                  case_when(strategy1=="Resource Acquisition" &
  #                              strategy2=="Resource Acquisition"~"Resource Acquisition",
  #                            strategy1=="Resource Use" &
  #                              strategy2=="Resource Use"~"Resource Use",
  #                            strategy1=="Stress Tolerance" &
  #                              strategy2=="Stress Tolerance"~"Stress Tolerance",
  #                            TRUE ~ "Between Strategies")) %>%
  #  dplyr::mutate(group = factor(`group`, levels = c("Resource Acquisition",
  #                                                   "Resource Use",
  #                                                   "Stress Tolerance",
  #                                                   "Between Strategies"), ordered = T))
  #p = ggplot2::ggplot(cor_melted1, aes(x=group, y=dist)) + geom_boxplot() +
  #  theme(axis.text=element_text(size=22),
  #        axis.title=element_blank())
  #ggplot2::ggsave(p, height = 6, width = 12,
  #       filename = file.path(outdir, paste0(dataset, ".alltraits.cov.boxplots.pdf")))

  #wilcox.res = wilcox.test(cor_melted1 %>% dplyr::filter(group == "within Resource Acquisition") %>% dplyr::pull(dist),
  #                         cor_melted1 %>% dplyr::filter(group == "between") %>% dplyr::pull(dist),
  #                         alternative = "greater")
  #wilcox.res = wilcox.test(cor_melted1 %>% dplyr::filter(group == "within Resource Use") %>% dplyr::pull(dist),
  #                         cor_melted1 %>% dplyr::filter(group == "between") %>% dplyr::pull(dist),
  #                         alternative = "greater")
  #wilcox.res = wilcox.test(cor_melted1 %>% dplyr::filter(group == "within Stress Tolerance") %>% dplyr::pull(dist),
  #                         cor_melted1 %>% dplyr::filter(group == "between") %>% dplyr::pull(dist),
  #                         alternative = "greater")
}

trait2trait_corr2 = function(feature_matrix) {
  temp = feature_matrix %>% dplyr::select(-id) %>% as.data.frame()
  ngenomes = nrow(temp)
  traits = colnames(temp)
  traits = setdiff(traits,
                   c(traits[which(apply(temp, 2, sum) == 0)],
                     traits[which(apply(temp, 2, sum) == ngenomes)]))

  cor_matrix = matrix(nrow = length(traits), ncol = length(traits))
  rownames(cor_matrix) = traits
  colnames(cor_matrix) = traits

  for(i in 1:length(traits)) {
    cat(i, "\n")
    for(j in 1:length(traits)) {
      if(is.double(temp[,traits[i]]) & is.double(temp[,traits[j]])) {
        test.result <- cor.test(x=temp[,traits[i]], y=temp[,traits[j]], method = 'spearman')
        rho = test.result$estimate
        p = test.result$p.value
        cor_matrix[traits[i], traits[j]] = rho
      }
      if(is.integer(temp[,traits[i]]) & is.integer(temp[,traits[j]])) {
        test.result <- cor.test(x=temp[,traits[i]], y=temp[,traits[j]], method = 'spearman')
        rho = test.result$estimate
        p = test.result$p.value
        cor_matrix[traits[i], traits[j]] = rho
      }
      # continuous and binary: glass rank biserial correlation coefficient
      if(is.double(temp[,traits[i]]) & is.integer(temp[,traits[j]])) {
        #i=1; j=49
        #test.result = glm(temp[,traits[j]] ~ temp[,traits[i]],family=binomial("logit"))
        #or = odds.ratio(test.result)[[1]][1]
        wilcoxonRG = wilcoxonRG(x = temp[,traits[i]], g = temp[,traits[j]])
        cor_matrix[traits[i], traits[j]] = wilcoxonRG
      }
      if(is.integer(temp[,traits[i]]) & is.double(temp[,traits[j]])) {
        #i=1; j=49
        #test.result = glm(temp[,traits[i]] ~ temp[,traits[j]],family=binomial("logit"))
        #or = odds.ratio(test.result)[[1]][1]
        wilcoxonRG = wilcoxonRG(x = temp[,traits[j]], g = temp[,traits[i]])
        cor_matrix[traits[i], traits[j]] = wilcoxonRG
      }
    }
  }

  cor_melted = melt_dist(cor_matrix)
  pdf(height = 6, width = 12,
      file.path(base, paste0(dataset, ".", "trait_matrixatgranularity3", "_traitcorrhist.pdf")))
  h = hist(cor_melted[,"dist"], breaks = seq(-1, 1, 0.1), xlab = "trait-trait correlation", ylab = "count", main = "", cex.lab = 1.5,xaxt= "n")
  axis(side=1, at = c(seq(-1, 0, 0.1), seq(0.1, 1, 0.1)), labels = c(seq(-1, 0, 0.1), seq(0.1, 1, 0.1)))
  dev.off()
  write.table(cor_melted,
              file = file.path(base, paste0(dataset, ".", "trait_matrixatgranularity3", "_traitcorr.xls")),
              row.names = F, col.names = T, sep = "\t", quote = F)

  col3 = colorRampPalette(c('red', 'white', 'blue'))
  col5 = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))

  pdf(height = 40,
      width = 40,
      file.path(base, paste0(dataset, ".alltraits.cov.ward-ordered.nolegend.pdf")))
  corrplot(cor_matrix,
           order = 'hclust', hclust.method = "ward",
           # c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid")
           #order = 'AOE', # angular order of the eigenvectors-Friendly 2002
           method = 'color',
           col = col3(16),
           col.lim = c(range(-1, range(cor_matrix)[2])),
           #addrect = 10, rect.lwd = 2,
           type = 'upper', diag = FALSE,
           cl.pos = "b", cl.cex = 3,
           tl.pos = "td", tl.cex = 1, tl.offset = 3, tl.col = "black",
           mar = c(2, 2, 2, 2)
           #type = 'full', tl.pos='n'
  )
  dev.off()
}

#' Calculate interguild variance
#'
#' @param genomeset_results
#'
#' @importFrom vegan adonis2
#' @importFrom ggplot2 ggplot
#' @export
calc_intergenome_variance = function(genomeset_results, 
                                     trait_granularity = "trait_matrixatgranularity3", 
                                     method = "wishart", 
                                     normby = "genome_length", 
                                     binarytype = "logical", 
                                     clustering_method = "ward.D",
                                     outdir, 
                                     ncores = 1) {
  ##############
  ### 
  ##############
  genomeset_results_norm = genomeset_results %>% 
    trait.normalize(normby = normby)
  genomeset_results_norm = genomeset_results_norm %>% 
    convert_traitdatatype(binarytype = binarytype)
  distance_matrix = genomeset_results_norm[[trait_granularity]] %>% 
    calc_mixeddist(method = method, binarytype = binarytype)

  hclust_genomes = hclust(distance_matrix, method = "ward.D")
  nguilds = seq(2, attr(distance_matrix, "Size")-1, 1)
  adonis_results_temp =
    parallel::mclapply(2:length(nguilds),
                       function(i) {
                          v = cutree(hclust_genomes, nguilds[i])
                          genome2cluster = data.frame(guild = factor(v))
                          rownames(genome2cluster) = names(v)
                          adonis_results = vegan::adonis2(distance_matrix ~ guild, data = genome2cluster, perm = 1)
                          result = data.frame(`number of guilds` = nguilds[i], 
                                              `percent variance` = adonis_results["guild", "R2"]*100,
                                              check.names = F)
                       },
                      mc.cores = ncores)
  adonis_results = do.call("rbind", adonis_results_temp)
  p = ggplot(adonis_results, aes(x=`number of guilds`, y=`percent variance`)) +
    geom_line() + geom_point(size=0.8) +
    #xlim(1, 500) +
    geom_hline(yintercept = 60, color = "red", size = 0.3) +
    geom_hline(yintercept = 70, color = "red", size = 0.3) +
    geom_hline(yintercept = 90, color = "red", size = 0.3) +
    scale_x_continuous(breaks = c(seq(0, attr(distance_matrix, "Size")-1, by = 5), attr(distance_matrix, "Size"))) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    ggtitle("Trait variance across genomes") +
    xlab("number of guilds") +
    ylab("% explained variance") +
    theme(text = element_text(size = 22),
          axis.text = element_text(size = 22))
  pdf_outfile = file.path(outdir, "traitvariance_acrossgenomes.pdf")
  ggsave(p, height = 8, width = 8,
         filename = pdf_outfile)
  png_outfile = file.path(outdir, "traitvariance_acrossgenomes.png")
  ggsave(p, height = 8, width = 8,
         filename = png_outfile)
  results = list(distance_matrix = distance_matrix, 
                 adonis_results = adonis_results, 
                 pdf_outfile = pdf_outfile, 
                 png_outfile = png_outfile)
  results
}



#' Genome by trait matrix heatmap.
#'
#' @param trait_matrix
#'
#' @import pheatmap
#' @export
#trait_matrix = trait_matrixatgranularity3_binary
#idcol = "id"
#annot_cols = c("mingentime", "ogt")
#granularity = 3
#clustering_distance_rows = genomeset_distances
#clustering_distance_cols = "binary"
#width = unit(40, "inches")
#height = unit(14.14*9, "inches")
#heatmap_width = width*0.8
#heatmap_height = height*0.8
#row_dend_width = width*0.05
#column_dend_height = height*0.02
#rightannotation_width = width*0.1
#topannotation_height = height*0.02
#bottomannotation_height = height*0.005
#pdf = TRUE; verbose = TRUE
#base_dir = "/Users/ukaraoz/Work/microtrait/code/github/microtrait-out"
#dataset = "soilgenomes"

cluster_traitmatrix = function(trait_matrix,
                               idcol = "id",
                               annot_cols = c("mingentime", "ogt"),
                               granularity = 3,
                               cutree_rows = NA,
                               cutree_cols = NA,
                               clustering_distance_rows,
                               clustering_distance_cols,
                               width = width, height = height,
                               heatmap_width, heatmap_height,
                               row_dend_width, column_dend_height,
                               rightannotation_width,
                               topannotation_height, bottomannotation_height,
                               outdir, dataset, pdf = TRUE, verbose = TRUE) {
  # data that goes into the heatmap
  toplot = trait_matrix %>% dplyr::select(-idcol) %>% dplyr::select(-annot_cols) %>% as.data.frame()
  rownames(toplot) = trait_matrix %>% dplyr::pull(idcol)
  colnames = colnames(toplot) %>% tibble::as.tibble() %>% dplyr::rename(`microtrait_trait-name` = value)

  annotation_colors = list(`microtrait_trait-strategy` = c("Resource Acquisition" = "red",
                                                           "Resource Use" = "blue",
                                                           "Stress Tolerance" = "darkgreen"),
                           `mingentime` = "yellow",
                           `ogt` = "brown")
  # any annotation data for rows, mingentime and ogt by default
  annotation_row = trait_matrix %>% dplyr::select(annot_cols) %>% as.data.frame()
  annot_range1 = annotation_row %>% dplyr::select(annot_cols[1]) %>% range()
  annot_range2 = annotation_row %>% dplyr::select(annot_cols[2]) %>% range()

  # any annotation data for columns
  annotation_col = colnames %>%
    dplyr::left_join(traits_listbygranularity[[granularity]], keep = TRUE, by = c("microtrait_trait-name" = "microtrait_trait-name"), suffix = c("", ".y")) %>%
    dplyr::select(c(`microtrait_trait-name`, `microtrait_trait-strategy`))
  ## add prevalence data
  prevalence = compute.prevalence(trait_matrix,type="trait_matrixatgranularity3")
  annotation_col = annotation_col %>%
    dplyr::left_join(prevalence, by = c("microtrait_trait-name" = "microtrait_trait-name")) %>%
    tibble::column_to_rownames(var = "microtrait_trait-name")

  # start building the pieces
  #GLOBAL_PADDING = unit(rep(0, 4), "points")
  ## top annotation for columns
  ht_opt("ROW_ANNO_PADDING" = unit(0.2, "inches"),
         "COLUMN_ANNO_PADDING" = unit(0.2, "inches"),
         heatmap_border = TRUE)
  p_top = ComplexHeatmap::HeatmapAnnotation(prevalence = anno_barplot(annotation_col[, "percent"],
                                                                      gp = gpar(col = NA, fill = "darkblue", lty = "blank"),
                                                                      axis_param = list(at = seq(0,100,20),
                                                                                        gp = gpar(fontsize = 16))),
                                            gp = gpar(col = "darkblue"),
                                            annotation_height = topannotation_height,
                                            annotation_label = "Trait Prevalence (%genomes positive)",
                                            annotation_name_gp = gpar(fontsize = 30),
                                            show_legend = FALSE)
  p_bottom = ComplexHeatmap::HeatmapAnnotation(strategy = annotation_col[, "microtrait_trait-strategy"],
                                               col = list(strategy = c("Resource Acquisition" = "red",
                                                                       "Resource Use" = "green",
                                                                       "Stress Tolerance" = "blue")),
                                               gp = gpar(col = "black"),
                                               annotation_height = bottomannotation_height,
                                               annotation_name_side = "left",
                                               annotation_name_gp = gpar(fontsize = 30),
                                               show_legend = FALSE)
  ## rigth annotation for rows
  p_right = ComplexHeatmap::rowAnnotation(
    mingentime = anno_barplot(annotation_row[, "mingentime"],
                              ylim = c(0, ceiling(annot_range1[2])),
                              gp = gpar(col = NA, fill = "darkblue", lty = "blank"),
                              axis_param = list(at = c(0, 2, 5, 10, 15, ceiling(annot_range1[2])),
                                                gp = gpar(fontsize = 16),
                                                labels = c("0", "2", "5", "10", "15", ceiling(annot_range1[2])),
                                                side = "bottom", facing = "outside",
                                                labels_rot = 45),
                              width = rightannotation_width/2),
    ogt = anno_points(annotation_row[, "ogt"],
                      ylim = c(0, ceiling(annot_range2[2])),
                      extend = 0.05,
                      gp = gpar(col = "darkblue", fontsize = 8),
                      axis_param = list(at = c(0, 20, 40, ceiling(annot_range2[2])),
                                        gp = gpar(fontsize = 16),
                                        labels = c("0", "20", "40", ceiling(annot_range2[2])),
                                        side = "bottom", facing = "outside",
                                        labels_rot = 45),
                      width = rightannotation_width/2),
    show_legend = FALSE,
    gap = unit(20, "points"),
    annotation_width = rightannotation_width,
    annotation_label = c("Minimal Doubling Time (hours)",
                         "Optimum Growth Temperature (C)"),
    annotation_name_gp = gpar(fontsize = 30),
    annotation_name_rot = 90)
  p_main = ComplexHeatmap::Heatmap(toplot, name = "main",
                                   clustering_distance_rows = clustering_distance_rows,
                                   clustering_distance_columns = clustering_distance_cols,
                                   col = c("white", "red"),
                                   use_raster = FALSE,
                                   show_row_names = F, show_column_names = T,
                                   row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 11),
                                   top_annotation = p_top,
                                   bottom_annotation = p_bottom,
                                   right_annotation = p_right,
                                   heatmap_height = heatmap_height,
                                   heatmap_width = heatmap_width,
                                   row_dend_width = row_dend_width,
                                   column_dend_height = column_dend_height,
                                   show_heatmap_legend = FALSE)
  #saveRDS(toplot, file.path(base, paste0(dataset, ".", "toplot.rds")))
  #pheatmap.outlist = list()
  # all: k=277; 70.19208
  # soil: k=194; 70.01773
  heatmap_outfile = file.path(base_dir, paste0(dataset, ".clustered_traitmatrix.pdf"))
  if(pdf == TRUE) {
    pdf(file = heatmap_outfile, width = width, height = height)
  }
  draw(p_main) #, padding = unit(c(1, 1, 1, 1), "points"))

  # add reference lines
  ComplexHeatmap::decorate_annotation("mingentime", {
    #grid.text("Age", unit(8, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
    grid.lines(x = unit(c(2, 2), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 5, col = "#7F000000"))
    grid.lines(x = unit(c(5, 5), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 5, col = "#7F000000"))
    grid.lines(x = unit(c(10, 10), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 5, col = "#7F000000"))
    grid.lines(x = unit(c(ceiling(annot_range1[2]), ceiling(annot_range1[2])), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 5, col = "#7F000000"))
  })

  ComplexHeatmap::decorate_annotation("ogt", {
    #grid.text("Age", unit(8, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
    grid.lines(x = unit(c(20, 20), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 5, col = "#7F000000"))
    grid.lines(x = unit(c(30, 30), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 5, col = "#7F000000"))
    grid.lines(x = unit(c(40, 40), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 5, col = "#7F000000"))
    grid.lines(x = unit(c(ceiling(annot_range2[2]), ceiling(annot_range2[2])), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 5, col = "#7F000000"))
  })

  ComplexHeatmap::decorate_annotation("prevalence", {
    #grid.text("Age", unit(8, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
    grid.lines(unit(c(0, 1), "npc"), unit(c(10, 10), "native"), gp = gpar(lty = 5, col = "#7F000000"))
    grid.lines(unit(c(0, 1), "npc"), unit(c(50, 50), "native"), gp = gpar(lty = 5, col = "#7F000000"))
    grid.lines(unit(c(0, 1), "npc"), unit(c(90, 90), "native"), gp = gpar(lty = 5, col = "#7F000000"))
  })
  if(pdf == TRUE) {
    dev.off()
  }

  # earlier attempt, pheatmap has many limitations
  #pheatmapout = pheatmap::pheatmap(toplot,
                                   #                       #kmeans_k = k, iter.max = 10, nstart = 25,
                                   #                       annotation_row = annotation_row,
                                   #                       annotation_col = annotation_col,
                                   #                       annotation_colors = annotation_colors,
                                   #                       cutree_rows = cutree_rows, cutree_cols = cutree_cols,
                                   #clustering_distance_rows = clustering_distance_rows,
                                   #clustering_distance_cols = clustering_distance_cols,
                                   #                       color = c("white", "red"), breaks = c(0,0.9,1),
                                   #                       show_rownames = F, show_colnames = T,
                                   #                       treeheight_row = 150, treeheight_col = 200,
                                   #                       #cellwidth = 1, cellheight = 1, labels_row = NA, border_color = NA, treeheight_row = 50,
                                   #                       height = height, width = width,
                                   #                       fontsize_row = 6, fontsize_col = 8, angle_col = 90, border_color = "black",
                                   #                       cluster_cols = T,
                                   #                       #annotation_row = annotation_row, annotation_colors = annotation_colors,
                                   #                       silent = T, legend = F,
                                   #                       filename = heatmap_outfile)
  #)
  if(verbose) {message("Clustered trait matrix plot: ", heatmap_outfile, ".")}
  #pheatmapout
}


#' Define guild matrix.
#'
#' @param trait_matrix trait_matrix
#' @param clusters_traitmatrix clusters_traitmatrix
#' @param nguilds nguilds
#'
#' @import pheatmap
#' @importFrom gtools stars.pval
#' @export
define_guilds = function(trait_matrix,
                         hclust_rows,
                         hclust_cols,
                         idcol = "id",
                         annot_cols = c("mingentime", "ogt"),
                         nguilds,
                         guildsizecutoff = 50,
                         traitvariability_threshold = 0.01,
                         width, height,
                         heatmap_width, heatmap_height,
                         rightannotation_width,
                         topannotation_height, bottomannotation_height,
                         outdir, dataset, verbose = TRUE) {
  #which.min(abs(maov_results[,2]-70))
  #hclust_rows = clusters_traitmatrix[[1]]
  hcut = cutree(hclust_rows, nguilds)
  genome2guild = data.frame(genome = names(hcut), guild = paste0("guild_", hcut)) %>% as_tibble()
  #  tibble::column_to_rownames(var = "microtrait_trait-name")

  temp = table(genome2guild[, "guild"])
  guild2size = data.frame(guild = names(temp), `numberofgenomes` = as.numeric(temp)) %>% as_tibble()
  genome2guild = genome2guild %>%
    dplyr::left_join(guild2size, by = c("guild" = "guild"))

  ############
  # plot guild size distribution
  outfile = file.path(outdir, paste0(dataset, ".guildsizedist.pdf"))
  temp = guild2size %>% pull(`numberofgenomes`) %>% table()
  guildsizedist = data.frame(guildsize = as.numeric(names(temp)), numberofguilds = as.numeric(temp))

  p = ggplot2::ggplot(guildsizedist,
                      aes(x=`guildsize`, y=`numberofguilds`)) +
      geom_bar(stat = "identity") +
      #scale_x_continuous(breaks = c(1,5,10,20,30,40,50, guildsizedist[which(guildsizedist[, "guildsize"] > 50), "guildsize"])) +
      scale_y_continuous(breaks = c(seq(0, max(guildsizedist[, "numberofguilds"]), 10), max(guildsizedist[, "numberofguilds"]))) +
      labs(x = "guild size (number of genomes)",
           y = "number of guilds") +
      theme(axis.text=element_text(size=12),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title=element_text(size=18,face="bold"))
  ggsave(p, height = 6, width = 12, filename = outfile)
  if(verbose) {message("Guild size distribution plot: ", outfile, ".")}
  ############

  ############
  # compute guild profiles
  count_traits = traits_listbygranularity[[as.numeric(granularities[3])]] %>%
    dplyr::filter(`microtrait_trait-type` == "count") %>%
    dplyr::pull(`microtrait_trait-name`) %>%
    as.character()

  binary_traits = traits_listbygranularity[[as.numeric(granularities[3])]] %>%
    dplyr::filter(`microtrait_trait-type` == "binary") %>%
    dplyr::filter(`microtrait_trait-name` != "Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation") %>%
    dplyr::pull(`microtrait_trait-name`) %>%
    as.character()

  guild2traitprofile = trait_matrix %>%
    dplyr::left_join(genome2guild, by = c("id" = "genome")) %>%
    dplyr::select(c("id", "guild", "numberofgenomes", everything())) %>%
    dplyr::select(-c("id", `numberofgenomes`)) %>%
    dplyr::group_by(guild) %>%
    dplyr::left_join(guild2size, by = c("guild" = "guild")) %>%
    dplyr::summarise(across(everything(), mean)) %>%
    dplyr::mutate_at(count_traits,
                     funs((. - min(.)) / (max(.) - min(.)))) %>%
    dplyr::arrange(desc(numberofgenomes)) %>%
    dplyr::select(c("guild", "numberofgenomes", everything())) %>% as.data.frame()
  # variability of positivity
  guild2variability = guild2traitprofile %>%
    dplyr::select(-c("guild", "numberofgenomes", annot_cols)) %>%
    dplyr::summarise(across(everything(), sd)) %>% t() %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "trait") %>%
    dplyr::rename(sd = V1)

  ############

  ############
  # build annotations for all guilds
  annotations = trait_matrix %>%
    dplyr::left_join(genome2guild, by = c("id" = "genome")) %>%
    dplyr::select(c("id", "guild", "numberofgenomes", everything())) %>%
    #dplyr::filter(numberofgenomes >= guildsizecutoff) %>%
    dplyr::arrange(desc(numberofgenomes)) %>%
    dplyr::select(id, guild, annot_cols)

  allguilds = guild2traitprofile %>% dplyr::pull(guild)
  mingentime_list = list(); ogt_list = list()
  for(i in 1:length(allguilds)) {
    mingentime_list[[allguilds[i]]] = annotations %>%
      dplyr::filter(guild == allguilds[[i]]) %>%
      dplyr::pull(`mingentime`)

    ogt_list[[allguilds[i]]] = annotations %>%
      dplyr::filter(guild == allguilds[[i]]) %>%
      dplyr::pull(`ogt`)
  }

  boxplots_allguilds = ComplexHeatmap::rowAnnotation(
    mingentime = anno_boxplot(mingentime_list,
                              ylim = c(0, mingentime_max),
                              #gp = gpar(col = NA, fill = "grey30", lty = "blank"),
                              axis_param = list(at = c(0, 2, 5, 10, 15, mingentime_max),
                                                labels = c("0", "2", "5", "10", "15", mingentime_max),
                                                side = "bottom", facing = "outside",
                                                labels_rot = 45),
                              width = rightannotation_width/2),
    ogt = anno_boxplot(ogt_list,
                       ylim = c(0, ogt_max),
                       #gp = gpar(col = NA, fill = "grey30", lty = "blank"),
                       axis_param = list(at = c(0, 20, 40, ogt_max),
                                         labels = c("0", "20", "40", ogt_max),
                                         side = "bottom", facing = "outside",
                                         labels_rot = 45),
                       width = rightannotation_width/2),
    show_legend = FALSE,
    gap = unit(5, "points"),
    annotation_label = c("Minimal Doubling Time (hours)",
                         "Optimum Growth Temperature (C)"),
    annotation_name_gp = gpar(fontsize = 10)
  )
  ############


  #guild2traitprofile_selectguilds = trait_matrix_bin %>%
  #  dplyr::left_join(genome2guild, by = c("id" = "genome")) %>%
  #  dplyr::select(c("id", "guild", "numberofgenomes", everything())) %>%
  #  dplyr::filter(numberofgenomes >= guildsizecutoff) %>%
  #  dplyr::arrange(desc(numberofgenomes))
#
  #guild2traitprofile_selectguilds_annotations = guild2traitprofile_selectguilds %>%
  #  dplyr::select(id, guild, annot_cols)

  # annotations
  ## pull trait names in the default order
  all_traits = traits_listbygranularity[[3]] %>%
      # dplyr::filter(`microtrait_trait-type` == "count") %>%
      dplyr::pull(`microtrait_trait-name`) %>%
      as.character()
  all_traits = all_traits %>% intersect(colnames(guild2traitprofile))

  selectguilds = guild2size %>%
    dplyr::arrange(desc(numberofgenomes)) %>%
    dplyr::filter(numberofgenomes >= guildsizecutoff) %>%
    dplyr::pull(guild)

  guild2traitprofile_selectguilds = guild2traitprofile %>%
    dplyr::filter(guild %in% selectguilds) %>%
    tibble::column_to_rownames(var = "guild")


  # significance testing
  trait_matrixforsign = trait_matrix %>%
    dplyr::left_join(genome2guild, by = c("id" = "genome")) %>%
    dplyr::select(c("id", "guild", "numberofgenomes", everything())) %>%
    dplyr::select(-c("id", `numberofgenomes`)) %>%
    dplyr::filter(guild %in% selectguilds)
  trait_matrix_p = matrix(nrow = 0, ncol = 2)
  for(c in 1:length(count_traits)) {
    temp = trait_matrixforsign %>%
      dplyr::select(c("guild", count_traits[c])) %>%
      dplyr::mutate(guild = factor(guild)) %>%
      as.data.frame
    pvalue = kruskal.test(get(colnames(temp)[2]) ~ guild, data = temp)$p.value
    effsize = temp %>% rstatix::kruskal_effsize(get(colnames(temp)[2]) ~ guild) %>% dplyr::select(effsize)
    trait_matrix_p = rbind(trait_matrix_p,
                           c(count_traits[c], pvalue, effsize))
    #cat(c, "\t", pvalue, "\n")
  }
  for(c in 1:length(binary_traits)) {
    temp = trait_matrixforsign %>%
      dplyr::select(c("guild", binary_traits[c])) %>%
      dplyr::mutate(guild = factor(guild)) %>%
      as.data.frame
    table = t(table(temp))
    pvalue = chisq.test(table)$p.value
    trait_matrix_p = rbind(trait_matrix_p,
                           c(binary_traits[c], pvalue))
    #cat(c, "\t", pvalue, "\n")
  }
  colnames(trait_matrix_p) = c("trait", "pvalue")
  trait_matrix_p = data.frame(trait_matrix_p,
                              pvalue.star = gtools::stars.pval(as.numeric(trait_matrix_p[,"pvalue"])))



  # right side annotations
  mingentime_max = ceiling(max(unlist(lapply(mingentime_list[selectguilds], max))))
  ogt_max = ceiling(max(unlist(lapply(ogt_list[selectguilds], max))))
  boxplots_selectguilds = ComplexHeatmap::rowAnnotation(
    mingentime = anno_boxplot(mingentime_list[selectguilds],
                              ylim = c(0, mingentime_max),
                              gp = gpar(col = "black", fill = "grey30"),
                              axis_param = list(at = c(0, 2, 5, 10, 15, mingentime_max),
                                                labels = c("0", "2", "5", "10", "15", mingentime_max),
                                                side = "bottom", facing = "outside",
                                                labels_rot = 45,
                                                gp=gpar(fontsize=16)),
                              width = rightannotation_width/2),
    ogt = anno_boxplot(ogt_list[selectguilds],
                       ylim = c(0, ogt_max),
                       gp = gpar(col = "black", fill = "grey30"),
                       axis_param = list(at = c(0, 20, 40, ogt_max),
                                         labels = c("0", "20", "40", ogt_max),
                                         side = "bottom", facing = "outside",
                                         labels_rot = 45,
                                         gp=gpar(fontsize=16)),
                       width = rightannotation_width/2),
    show_legend = FALSE,
    gap = unit(5, "points"),
    annotation_label = c("Minimal Doubling Time (hours)",
                         "Optimum Growth Temperature (C)"),
    annotation_name_gp = gpar(fontsize = 20)
  )

  toplot = guild2traitprofile_selectguilds
  selecttraits = guild2variability %>% dplyr::filter(sd > traitvariability_threshold) %>% dplyr::pull(`trait`)
  toplot = toplot %>%
    dplyr::select(selecttraits) %>% as.data.frame()
  row_labels = paste0(rownames(toplot), " (",
                      rownames(toplot) %>% as_tibble() %>%
                        dplyr::rename(guild = value) %>%
                        dplyr::left_join(guild2size) %>%
                        dplyr::pull(numberofgenomes), " genomes)")

  p_top = ComplexHeatmap::HeatmapAnnotation(`variability` = anno_barplot(guild2variability[guild2variability$trait %in% selecttraits, "sd"],
                                                                         gp = gpar(col = "black", fill = "darkblue", lty = "blank"),
                                                                         axis_param = list(at = seq(0,0.5,0.1),
                                                                                           gp = gpar(fontsize = 15))
                                                                         ),
                                            gp = gpar(col = "black"),
                                            annotation_height = topannotation_height,
                                            annotation_label = "variability (stdev)",
                                            annotation_name_gp = gpar(fontsize = 16),
                                            show_legend = FALSE)
  p_top_sign = ComplexHeatmap::HeatmapAnnotation(
    pvalue = trait_matrix_p[trait_matrix_p$trait %in% selecttraits, "pvalue.star"],
          #`pvalue` = anno_text(trait_matrix_p[trait_matrix_p$trait %in% selecttraits, "pvalue.star"],
          #                                                            rot = 90,
          #                                                            which = "column",
          #                                                            gp = gpar(fontsize = 9)
          #                                                           ),
    col = list(pvalue = c("*" = "#FEE0D2",
                            "**" = "#EF3B2C",
                            "***" = "#67000D",
                            " " = "gray")),
    gp = gpar(col = "black"),
    annotation_height = topannotation_height,
    annotation_label = "pvalue",
    annotation_name_gp = gpar(fontsize = 16),
    show_legend = TRUE)

  p_bottom = ComplexHeatmap::HeatmapAnnotation(
    strategy = colnames(toplot) %>% as_tibble() %>%
                dplyr::rename(`microtrait_trait-name` = value) %>%
                dplyr::left_join(traits_listbygranularity[[3]], keep = TRUE, by = c("microtrait_trait-name" = "microtrait_trait-name"), suffix = c("", ".y")) %>%
                dplyr::pull(`microtrait_trait-strategy`),
    col = list(strategy = c("Resource Acquisition" = "red",
                            "Resource Use" = "green",
                            "Stress Tolerance" = "blue")),
    gp = gpar(col = "black"),
    annotation_height = bottomannotation_height,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 16),
    show_legend = FALSE)


  p_main = ComplexHeatmap::Heatmap(toplot,
                                   col = colorRamp2(c(0, 1), c("white", "#CB181D")),
                                   #col = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(20),
                                   cluster_rows = FALSE, cluster_columns = FALSE,
                                   show_row_names = TRUE, show_column_names = TRUE,
                                   show_row_dend = FALSE, show_column_dend = FALSE,
                                   show_heatmap_legend = TRUE,
                                   row_labels = row_labels,
                                   row_names_side = "left",
                                   row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 14),
                                   heatmap_height = heatmap_height,
                                   heatmap_width = heatmap_width,
                                   #row_split = rownames(guild2traitprofile_selectguilds_toplot[, trait_order]),
                                   top_annotation = p_top_sign,
                                   bottom_annotation = p_bottom,
                                   right_annotation = boxplots_selectguilds
                                   )

  heatmap_outfile = file.path(base_dir, paste0(dataset, ".guild2traitprofilewlegend.pdf"))
  if(pdf == TRUE) {
    pdf(file = heatmap_outfile, width = width, height = height)
  }
  ht_opt("ROW_ANNO_PADDING" = unit(0.05, "inches"),
         "COLUMN_ANNO_PADDING" = unit(0.05, "inches"),
         heatmap_border = TRUE)
  draw(p_main, padding = unit(c(10, 0, 0, 0), "inches"))

  # add reference lines
  ComplexHeatmap::decorate_annotation("mingentime", {
    #grid.text("Age", unit(8, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
    grid.lines(x = unit(c(2, 2), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
    grid.lines(x = unit(c(5, 5), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
    grid.lines(x = unit(c(10, 10), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
    grid.lines(x = unit(c(mingentime_max, mingentime_max), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
  })
  ComplexHeatmap::decorate_annotation("ogt", {
    #grid.text("Age", unit(8, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
    grid.lines(x = unit(c(20, 20), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
    grid.lines(x = unit(c(30, 30), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
    grid.lines(x = unit(c(40, 40), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
    grid.lines(x = unit(c(ogt_max, ogt_max), "native"), y = unit(c(0, 1), "npc"), gp = gpar(lty = 2, col = "black"))
  })
  dev.off()


  # outfile = file.path(outdir, paste0(dataset, ".selectguild.profiles.pdf"))
  # pheatmap::pheatmap(guild2traitprofile_selectguilds_toplot[, trait_order],
  #          color=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Reds"))(20),
  #          width = 44, height = 22,
  #          cellwidth = 16, cellheight = 16,
  #          fontsize_row = 6, fontsize_col = 14, angle_col = 90, border_color = "black",
  #          cluster_rows = FALSE, cluster_cols = FALSE,
  #          gaps_row = 1:nrow(trait_matrix_bin_selectguilds_toplot),
  #          legend_breaks = seq(0,1,0.1),
  #          legend_labels = paste0("trait positivity: ", as.character(seq(0,1,0.1))),
  #          show_rownames = T, show_colnames = T,
  #          filename = outfile)
  if(verbose) {message("Guild profile plot for guilds larger than ", guildsizecutoff,
                       " genomes: ", heatmap_outfile, ".")}
  result = list()
  result[["genome2guild"]] = genome2guild
  result[["guild2traitprofile"]] = guild2traitprofile
  guild2traitprofile_outfile = file.path(base_dir, paste0(dataset, ".guild2traitprofile.xls"))
  write.table(result[["guild2traitprofile"]][, c("guild", "numberofgenomes", selecttraits)],
              file = guild2traitprofile_outfile,
              row.names = F, col.names = T, sep = "\t", quote = F)
  if(verbose) {message("Guild profiles for all guilds: ", guild2traitprofile_outfile, ".")}

  result
}

#' Melt a square distance matrix into long format
#'
#' This will take a square distance matrix, and will transform in to long
#' format. It will remove upper triangle, and diagonal elements, so you
#' end with only (n)*(n-1)/2 rows, where n are the total number of rows in
#' the distance matrix.
#'
#' @param dist An object of class matrix, it must be square
#' @param order A character vector of size n with the order of the columns and rows (default: NULL)
#' @param dist_name A string to name the distance column in the output (default: dist)
#'
#' @return A data.frame with three columns: (1) iso1; (2) iso2; (3) dist. iso1 and
#' iso2 indicate the pair being compared, and dist indicates the distance between
#' that pair.
#'
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom lazyeval interp
#' @export
melt_dist <- function(dist, order = NULL, dist_name = 'dist') {
  if(!is.null(order)){
    dist <- dist[order, order]
  } else {
    order <- row.names(dist)
  }
  diag(dist) <- NA
  dist[upper.tri(dist)] <- NA
  dist_df <- as.data.frame(dist)
  dist_df$variable1 <- row.names(dist)
  dist_df <- dist_df %>%
    tidyr::gather_(key = "variable2", value = lazyeval::interp("dist_name", dist_name = as.name(dist_name)), order, na.rm = T)
  return(dist_df)
}

#' Count prevalence of traits, rules, or hmms
#'
#' @param genomes_matrix genomes_matrix
#' @param type type
#' @return results
#'
#' @export compute.prevalence
compute.prevalence <- function(feature_matrix, type, idcol = "id") {
  if(type == "hmm_matrix") {
    prevalence = feature_matrix %>% tibble::as_tibble() %>%
      dplyr::select(c(idcol, intersect(colnames(feature_matrix), hmms_fromrules %>% pull(`microtrait_hmm-name`)))) %>%
      tidyr::pivot_longer(!id, names_to = "hmm", values_to = "presence") %>%
      dplyr::select(-id) %>%
      dplyr::count(hmm, presence,.drop=FALSE) %>%
      dplyr::group_by(hmm) %>%
      dplyr::mutate(percent = n/sum(n)*100) %>%
      dplyr::filter(`presence` == 0) %>%
      dplyr::mutate(`percent` = 100-`percent`) %>%
      dplyr::select(c("hmm", "percent")) %>%
      dplyr::inner_join(hmms_fromrules, by = c("hmm" = "microtrait_hmm-name"))
  }

  if(type == "rule_matrix") {
    prevalence = feature_matrix %>% tibble::as_tibble() %>%
      dplyr::select(c(idcol, rules %>% pull(`microtrait_rule-name`))) %>%
      tidyr::pivot_longer(!id, names_to = "rule", values_to = "presence") %>%
      dplyr::select(-id) %>%
      dplyr::count(rule, presence,.drop=FALSE) %>%
      dplyr::group_by(rule) %>%
      dplyr::mutate(percent = n/sum(n)*100) %>%
      dplyr::filter(`presence` == 0) %>%
      dplyr::mutate(`percent` = 100-`percent`) %>%
      dplyr::select(c("rule", "percent")) %>%
      dplyr::inner_join(rules, by = c("rule" = "microtrait_rule-name"))
  }

  if(type == "trait_matrixatgranularity3") {
    prevalence = trait_matrix %>% tibble::as_tibble() %>%
      dplyr::select(c(idcol,
                      intersect(colnames(trait_matrix),
                                traits_listbygranularity[[3]] %>%
                                  #dplyr::filter(`microtrait_trait-type` == "binary") %>%
                                  dplyr::pull(`microtrait_trait-name`)))) %>%
      # conversion to factor required for all 0/all 1 traits
      dplyr::mutate_at(vars(!starts_with("id")),
                       funs(case_when(. >= 1 ~ factor(1, levels = c(0, 1)),
                                      TRUE ~ factor(0, levels = c(0, 1))))) %>%
      tidyr::pivot_longer(!idcol, names_to = "microtrait_trait-name", values_to = "presence") %>%
      dplyr::select(-idcol) %>%
      dplyr::count(`microtrait_trait-name`, presence,.drop=FALSE) %>%
      dplyr::group_by(`microtrait_trait-name`) %>%  #methaneoxidation_genes = c("mmoX", "mmoY", "mmoZ", "mmoC", "mmoD", "pmoA-amoA", "pmoB-amoB", "pmoC-amoC")
      dplyr::mutate(percent = n/sum(n)*100) %>%
      dplyr::filter(`presence` == 1) %>%
      # dplyr::mutate(`percent` = 100-`percent`) %>%
      dplyr::select(c("microtrait_trait-name", "n", "percent"))
  }
  prevalence
}
