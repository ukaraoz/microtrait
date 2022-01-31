#' Build lme model
#'
#' @importFrom ape varcomp
#' @importFrom nlme lme
#' @return model
#'
#' @export build.taxa.model
build.taxa.model <- function(trait_matrix, trait, species = TRUE, optimizer = "nlm") {
  # spaces in variable names!
  trait_matrix = trait_matrix %>%
    dplyr::rename_with(~ gsub(" |:|-|,|/|\\(|\\)", "_", .x), all_of(trait))
  trait = gsub(" |:|-|,|/|\\(|\\)", "_", trait)

  if(species == TRUE) {
    model1.trait <- try(nlme::lme(as.formula(paste0(trait,"~1")),
                              random = ~1|NCBI_Phylum/NCBI_Class/NCBI_Order/NCBI_Family/NCBI_Genus/NCBI_Species/Ecosystem_Category,
                              data = trait_matrix,
                              na.action = na.omit,
                              control = nlmeControl(maxIter = 200, msMaxIter = 200, opt=optimizer)),
                        silent = TRUE)
  } else {
    model1.trait <- try(nlme::lme(as.formula(paste0(trait,"~1")),
                                  random = ~1|NCBI_Phylum/NCBI_Class/NCBI_Order/NCBI_Family/NCBI_Genus/Ecosystem_Category,
                                  data = trait_matrix,
                                  na.action = na.omit,
                                  control = nlmeControl(maxIter = 200, msMaxIter = 200, opt=optimizer)),
                        silent = TRUE)
  }
  if(class(model1.trait) != "try-error") {
    model1.trait.var.temp = round(ape::varcomp(model1.trait,T,F)*100,5)
    model1.trait.var = data.frame(variable = names(model1.trait.var.temp[1:8]),
                                  variance = model1.trait.var.temp[1:8],
                                  trait = trait,
                                  model = paste0(as.character(model1.trait$call)[4])) %>% tibble::as_tibble()
  } else {
    model1.trait.var = data.frame(variable = rep("error", 7),
                                  variance = rep(0, 7),
                                  trait = trait,
                                  model = rep(model1.trait[1], 7)) %>% tibble::as_tibble()
  }


  if(species == TRUE) {
    model2.trait <- try(nlme::lme(as.formula(paste0(trait,"~1")),
                                  random = ~1|Ecosystem_Category/NCBI_Phylum/NCBI_Class/NCBI_Order/NCBI_Family/NCBI_Genus/NCBI_Species,
                                  data = trait_matrix,
                                  na.action = na.omit,
                                  control = nlmeControl(maxIter = 200, msMaxIter = 200, opt=optimizer)),
                        silent = TRUE)
  } else {
    model2.trait <- try(nlme::lme(as.formula(paste0(trait,"~1")),
                                  random = ~1|Ecosystem_Category/NCBI_Phylum/NCBI_Class/NCBI_Order/NCBI_Family/NCBI_Genus,
                                  data = trait_matrix,
                                  na.action = na.omit,
                                  control = nlmeControl(maxIter = 200, msMaxIter = 200, opt=optimizer)),
                        silent = TRUE)
  }

  if(class(model2.trait) != "try-error") {
    model2.trait.var.temp = try(round(varcomp(model2.trait,T,F)*100,5),
                                silent = TRUE)
    model2.trait.var = data.frame(variable = names(model2.trait.var.temp[1:8]),
                                  variance = model2.trait.var.temp[1:8],
                                  trait = trait,
                                  model = paste0(as.character(model2.trait$call)[4])) %>% tibble::as_tibble()
  } else {
    model2.trait.var = data.frame(variable = rep("error", 7),
                                  variance = rep(0, 7),
                                  trait = trait,
                                  model = rep(model2.trait[1], 7)) %>% tibble::as_tibble()
  }
  result = dplyr::bind_rows(model1.trait.var,
                            model2.trait.var)
  result
}


#' Run consentrait script
#'
#' @return result
#'
#' @export run.consentrait.singletrait
run.consentrait.singletrait <- function(trait, traitmatrix, tree, treeroot,
                                        nbootstrap = 100, nproc = 40, perc.share.cutoff = 90,
                                        outdir, consentraitbin = "/global/homes/u/ukaraoz/bin/consentrait.R",
                                        run = FALSE) {
  # traitmatrix here is for a single trait: ncol = 2
  tic("consentrait")
  trait_nospace = gsub(" |:|-|,|/|\\(|\\)", "__", trait)
  traitfile = file.path(outdir, paste0("traitmatrix_", trait_nospace, ".txt"))
  treefile = file.path(outdir, paste0("traitmatrix_", trait_nospace, ".tree"))
  clustersizefileprefix = paste0("clustersize_", trait_nospace)
  clusterdistfileprefix = paste0("clusterdist_", trait_nospace)
  taudtablefile = file.path(outdir, paste0("taud_", trait_nospace, ".txt"))
  taudbootstrapfile = file.path(outdir, paste0("taudbootstrap_", trait_nospace, ".txt"))
  outfile = file.path(outdir, paste0("consentrait.", trait_nospace, ".txt"))
  errfile = file.path(outdir, paste0("consentrait.", trait_nospace, ".err.txt"))

  write.table(trait_matrix.sample[, c("treetip", trait)],
              file = traitfile,
              sep = "\t",
              quote=FALSE,
              row.names = F,
              col.names = F)
  write.tree(tree.sample,
             file = treefile)

  command = paste0("Rscript", " ", consentraitbin, " ",
                   "-b", " ", nbootstrap, " ",
                   "-p", " ", nproc, " ",
                   "-s", " ", perc.share.cutoff, " ",
                   "-c", " ", clustersizefileprefix, " ",
                   "-d", " ", clusterdistfileprefix, " ",
                   "-t", " ", taudtablefile, " ",
                   "-u", " ", taudbootstrapfile, " ",
                   treefile, " ",
                   treeroot, " ",
                   traitfile,
                   " > ", outfile, " > ", errfile)


  if(run == TRUE) {
    message("Executing: ", command, "\n")
    system(command)
    toc(log = TRUE, quiet = TRUE)
    log.txt <- tic.log(format = TRUE)
    result = data.frame(command = command, time = log.txt[[1]])
    return(result)
  } else {
    message("Executing: ", command, "\n")
  }
}
