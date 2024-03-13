library(microtrait)
library(dplyr)

# Essentially, a list of commands on the main markdown document
# collect isolate results
#base = "/global/homes/u/ukaraoz/cscratch/alltarballs"
base_dir = "/Users/ukaraoz/Work/microtrait/code/github/microtrait-out"
dataset = "environmentalgenomes"
#rds_files = list.files(file.path(base_dir, "microtrait-out"),
#                       full.names = T, recursive = F, pattern = ".microtrait.rds$")
# genomeset_results = make.genomeset.results(rds_files = rds_files,
#                                            ids = sub(".microtrait.rds", "", basename(rds_files)),
#                                            ncores = floor(0.8*detectCores()))
# saveRDS(genomeset_results, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
genomeset_results = readRDS(file.path(base_dir, paste0(dataset, ".microtraitresults.backup.rds")))
# add metadata and gp
genome_metadata = readRDS(file.path(base_dir, paste0(dataset, ".metadata.rds")))
#write.table(genome_metadata, file.path(base_dir, "environmentalgenomes.metadata.xls"), row.names=F, col.names=T, sep ="\t", quote =F)
genomeset_results_wmetadata = add.metadata(genomeset_results, genome_metadata, genome_metadata_idcol = "IMG Taxon ID")

gp_results = readRDS(file.path(base_dir, paste0(dataset, ".gp.rds"))) %>%
  dplyr::select(-c("sdgentime", "nNonHEG", "nHEG", "OGT")) # %>%
  #dplyr::rename(nHEG.VieiraSilva = nHEG) %>%
  #dplyr::rename(mingentime.VieiraSilva = mingentime) %>%
  #dplyr::rename(ogt.Zeldovich = OGT)

#genomeset_results_wmetadata = add.metadata(genomeset_results_wmetadata, gp_results, genome_metadata_idcol = "id")

grodon.results = combine.grodon.results(list.files("/Users/ukaraoz/Work/microtrait/code/inst/extdata/grodon", pattern = "_growth.rds$", full.names = T)) %>%
  dplyr::select(c("id", "mingentime")) # %>%
  #dplyr::rename(nHEG = nhighlyexpressed)
#genomeset_results_wmetadata = add.metadata(genomeset_results_wmetadata, grodon.results, genome_metadata_idcol = "id")

ogt.sauer.results = read.table("/Users/ukaraoz/Work/microtrait/code/inst/extdata/organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID_taxonids.OGT.txt", header = T, sep = "\t") %>%
  as.tibble() %>%
  dplyr::select(c("id", "ogt")) %>%
  dplyr::mutate(id = as.character(id)) %>%
  dplyr::rename(ogt = ogt)

for(i in 1:3) {
  genomeset_results[[i]] = genomeset_results[[i]] %>%
    dplyr::left_join(grodon.results, by = c("id" = "id")) %>%
    dplyr::left_join(ogt.sauer.results, by = c("id" = "id"))
}
saveRDS(genomeset_results, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))

genomeset_results_wmetadata = add.metadata(genomeset_results_wmetadata, ogt.sauer.results, genome_metadata_idcol = "id")
genomeset_results_wmetadata = genomeset_results_wmetadata %>% convert_traitdatatype(binarytype = "logical")
#saveRDS(genomeset_results_wmetadata, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
genomeset_results_wmetadata = readRDS(file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))

# grodon.results = combine.grodon.results(list.files("/Users/ukaraoz/Work/microtrait/code/inst/extdata/grodon", pattern = "_growth.rds$", full.names = T))
# grodon.results = grodon.results %>%
#   dplyr::filter(mingentime != "NA") %>%
#   dplyr::filter(nhighlyexpressed > 20) %>%
#   dplyr::filter(nhighlyexpressed < 50)
# gp_results = gp_results %>%
#   dplyr::filter(nHEG > 20) %>%
#   dplyr::filter(nHEG < 50)
# temp = grodon.results %>%
#   dplyr::inner_join(gp_results, by = c("id" = "id"))
# a = temp %>% pull("mingentime.x")
# b = temp %>% pull("mingentime.y")
##################
# normalize by genome size and code binary traits as factor
genomeset_results_wmetadata_norm =
  genomeset_results_wmetadata %>%
  microtrait::trait.normalize(normby = "Estimated Size")

# filter to soil genomes, i.e. exclude aquatic
traits = traits_listbygranularity[[3]] %>%
  dplyr::select(`microtrait_trait-name`) %>%
  dplyr::filter(`microtrait_trait-name` != "Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation") %>%
  dplyr::pull(`microtrait_trait-name`)
trait_matrixatgranularity3 = genomeset_results_wmetadata_norm[["trait_matrixatgranularity3"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  #dplyr::select(-`Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation`) %>%
  dplyr::select(c("id", traits, "mingentime", "ogt")) %>%
  #dplyr::slice(1:100) %>%
  dplyr::filter(`mingentime` < 100) %>%   # max out mingentime at 25 days
  dplyr::filter(!is.na(`mingentime`) & !is.na(`ogt`))
trait_matrixatgranularity3_binary = trait_matrixatgranularity3 %>% microtrait::trait.continuous2binary()

trait2trait_corr(trait_matrixatgranularity3_binary, verbose = TRUE, idcol = "id", outdir = base_dir, dataset = "soilgenomes")

####
genomeset_distances = trait_matrixatgranularity3 %>% microtrait::calc_mixeddist(idcol = "id", col2ignore = c("mingentime", "ogt"), method = "wishart", binarytype = "logical", byrow = 1, verbose = TRUE)
prevalence = compute.prevalence(trait_matrixatgranularity3_binary, type="trait_matrixatgranularity3")

A4_ratio = 11.75/8.25
width = unit(8.25*4, "inches")
height = unit(8.25*4*A4_ratio, "inches")
cluster_traitmatrix(trait_matrix = trait_matrixatgranularity3_binary,
                    idcol = "id", annot_cols = c("mingentime", "ogt"), granularity = 3,
                    clustering_distance_rows = genomeset_distances, clustering_distance_cols = "binary",
                    width = width, height = height,
                    heatmap_width = width*0.8, heatmap_height = height*0.70,
                    row_dend_width = width*0.05, column_dend_height = height*0.02,
                    rightannotation_width = width*0.1, topannotation_height = height*0.02,
                    bottomannotation_height = height*0.005,
                    outdir = base_dir, dataset = "soilgenomes", pdf = TRUE)

#load precomputed
maov_results_list = readRDS("/Users/ukaraoz/Work/microtrait/code/github/microtrait-out/soilgenomes_maov_results_list.rds")
maov_results = matrix(nrow = 0, ncol = 4)
nguilds = seq(2, 3556, 2)
for(i in 1:length(maov_results_list)) {
  #cat(i, "\n")
  maov_results = rbind(maov_results,
                       c(nguilds[i], maov_results_list[[i]]$R2))
}

maov_results = data.frame(`number of guilds` = maov_results[,1],
                          `percent variance` = round(maov_results[,2]*100,2),
                          check.names = F)
p = ggplot2::ggplot(maov_results,aes(x=`number of guilds`, y=`percent variance`)) +
  geom_line() + geom_point(size=0.8) +
  #xlim(1, 500) +
  geom_hline(yintercept = 60, color = "red", size = 0.3) +
  geom_vline(xintercept = 100, color = "red", size = 0.3) +
  geom_hline(yintercept = 70, color = "red", size = 0.3) +
  geom_vline(xintercept = 194, color = "red", size = 0.3) +
  geom_hline(yintercept = 80.04, color = "red", size = 0.3) +
  geom_vline(xintercept = 400, color = "red", size = 0.3) +
  scale_x_continuous(breaks = c(seq(0, 3556, by = 200), 3556)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  theme(axis.text=element_text(size=20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title=element_text(size=22,face="bold"))
dataset = "soilgenomes"
ggsave(p, height = 6, width = 12,
       filename = file.path(base_dir, paste0(dataset, ".guilds.percentvariance.pdf")))

nguildsat70perc = maov_results[which.min(abs(maov_results[,2]-70)), "number of guilds"]

toplot = trait_matrixatgranularity3_binary %>% dplyr::select(-idcol) %>% dplyr::select(-annot_cols) %>% as.data.frame()
rownames(toplot) = trait_matrix %>% dplyr::pull(idcol)
colnames = colnames(toplot) %>% as.tibble() %>% dplyr::rename(`microtrait_trait-name` = value)
pheatmapout = pheatmap::pheatmap(toplot,
                                 clustering_distance_rows = clustering_distance_rows,
                                 clustering_distance_cols = clustering_distance_cols)
hclust_rows = pheatmapout[[1]]
hclust_cols = pheatmapout[[2]]
nguilds = nguildsat70perc
guildsizecutoff = 50
width = unit(40, "inches"); height = unit(20, "inches")
heatmap_width = width*0.96; heatmap_height = height*0.45
topannotation_height = unit(1, "inches"); bottomannotation_height = unit(0.3, "inches")
rightannotation_width = width*0.2
traitvariability_threshold = 0.01
defined_guilds = microtrait::define_guilds(trait_matrixatgranularity3,
                                           hclust_rows = hclust_rows,
                                           hclust_cols = hclust_cols,
                                           nguilds = nguildsat70perc,
                                           guildsizecutoff = 50, traitvariability_threshold = 0.01,
                                           width = width, height = height,
                                           heatmap_width = heatmap_width, heatmap_height = heatmap_height,
                                           rightannotation_width = rightannotation_width,
                                           topannotation_height = topannotation_height,
                                           bottomannotation_height = bottomannotation_height,
                                           outdir = base_dir, dataset = "soilgenomes", verbose = TRUE)


library(microtrait)
library(dplyr)

# collect isolate results
#base = "/global/homes/u/ukaraoz/cscratch/alltarballs"
base_dir = "/Users/ukaraoz/Work/microtrait/code/github/test/microtrait-out"
dataset = "terrestrialgenomes"
rds_files = list.files(file.path(base_dir, "microtrait-out"),
                       full.names = T, recursive = F, pattern = ".microtrait.rds$")
genomeset_results = make.genomeset.results(rds_files = rds_files,
                                           ids = sub(".microtrait.rds", "", basename(rds_files)),
                                           ncores = floor(0.8*detectCores()))
saveRDS(genomeset_results, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))

# add metadata and gp
genome_metadata = readRDS(file.path(base_dir, paste0(dataset, ".metadata.rds")))
#write.table(genome_metadata, file.path(base, "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID.metadata.xls"), row.names=F, col.names=T, sep ="\t", quote =F)
genomeset_results_wmetadata = add.metadata(genomeset_results, genome_metadata, genome_metadata_idcol = "IMG Taxon ID")

gp_results = readRDS(file.path(base_dir, paste0(dataset, ".gp.rds"))) %>%
  dplyr::select(-c("sdgentime", "nHEG", "nNonHEG")) %>%
  dplyr::rename(mingentime.VieiraSilva = mingentime) %>%
  dplyr::rename(ogt.Zeldovich = OGT)
genomeset_results_wmetadata = add.metadata(genomeset_results_wmetadata, gp_results, genome_metadata_idcol = "id")

grodon.results = combine.grodon.results(list.files("/Users/ukaraoz/Work/microtrait/code/inst/extdata/grodon", pattern = "_growth.rds$", full.names = T)) %>%
  dplyr::select(c("id", "mingentime"))
genomeset_results_wmetadata = add.metadata(genomeset_results_wmetadata, grodon.results, genome_metadata_idcol = "id")

ogt.sauer.results = read.table("/Users/ukaraoz/Work/microtrait/code/inst/extdata/organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID_taxonids.OGT.txt", header = T, sep = "\t") %>%
  as.tibble() %>%
  dplyr::select(c("id", "ogt")) %>%
  dplyr::mutate(id = as.character(id)) %>%
  dplyr::rename(ogt = ogt)
genomeset_results_wmetadata = add.metadata(genomeset_results_wmetadata, ogt.sauer.results, genome_metadata_idcol = "id")
saveRDS(genomeset_results_wmetadata, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))

# grodon.results = combine.grodon.results(list.files("/Users/ukaraoz/Work/microtrait/code/inst/extdata/grodon", pattern = "_growth.rds$", full.names = T))
# grodon.results = grodon.results %>%
#   dplyr::filter(mingentime != "NA") %>%
#   dplyr::filter(nhighlyexpressed > 20) %>%
#   dplyr::filter(nhighlyexpressed < 50)
# gp_results = gp_results %>%
#   dplyr::filter(nHEG > 20) %>%
#   dplyr::filter(nHEG < 50)
# temp = grodon.results %>%
#   dplyr::inner_join(gp_results, by = c("id" = "id"))
# a = temp %>% pull("mingentime.x")
# b = temp %>% pull("mingentime.y")

##################

rds_dir = "/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download/PRJNA348753_Parks/microtrait"
out_dir = "/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download/PRJNA348753_Parks"
dataset = "PRJNA348753_Parks"


genomeset_results = readRDS(file.path(out_dir, paste0(dataset, ".microtraitresults.rds")))
write.genomeset.results(genomeset_results, out_dir, dataset)
ids = sub(".microtrait.rds", "", basename(rds_files))

rule_matrix3 = combine.results(rds_files[8001:length(rds_files)], type = "rule", ncores = 20) %>% tibble::as_tibble()
rule_matrix = bind_rows(rule_matrix1, rule_matrix2, rule_matrix3)


saveRDS(trait_matrixatgranularity1, file.path(out_dir, paste0(dataset, ".trait_matrixatgranularity1.rds")))
saveRDS(trait_matrixatgranularity2, file.path(out_dir, paste0(dataset, ".trait_matrixatgranularity2.rds")))
saveRDS(trait_matrixatgranularity3, file.path(out_dir, paste0(dataset, ".trait_matrixatgranularity3.rds")))
trait_matrixatgranularity1 = readRDS(file.path(out_dir, paste0(dataset, ".trait_matrixatgranularity1.rds")))
trait_matrixatgranularity2 = readRDS(file.path(out_dir, paste0(dataset, ".trait_matrixatgranularity2.rds")))
trait_matrixatgranularity3 = readRDS(file.path(out_dir, paste0(dataset, ".trait_matrixatgranularity3.rds")))
hmm_matrix = readRDS(file.path(out_dir, paste0(dataset, ".hmm_matrix.rds")))
rule_matrix = readRDS(file.path(out_dir, paste0(dataset, ".rule_matrix.rds")))

saveRDS(hmm_matrix, file.path(out_dir, paste0(dataset, ".hmm_matrix.rds")))
saveRDS(rule_matrix, file.path(out_dir, paste0(dataset, ".rule_matrix.rds")))

saveRDS(genomeset_results, file.path(out_dir, paste0(dataset, ".microtraitresults.rds")))

library(microtrait)
library(dplyr)
rds_dir = "/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download/PRJNA348753_Parks"
genomeset_results = make.genomeset.results(rds_dir, ncores)
saveRDS(genomeset_results, file.path("/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download", paste0("PRJNA348753_Parks", ".microtraitresults.rds")))

rds_dir = "/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download/PRJNA348753_Parks"
rds_files = list.files(rds_dir, full.names = T, recursive = F, pattern = ".*rds$")


###########

base = "/Users/ukaraoz/Work/microtrait/code/microtrait-analysis/img.rds/growthpred"
files = list.files(base, full.names = T, recursive = F, pattern = ".*growthpred.txt.results$")
ids = sub(".growthpred.txt.results", "", basename(files))
gp.results = combine.gp.results(files, ids)
saveRDS(gp.results, file.path("/Users/ukaraoz/Work/microtrait/code/microtrait/inst/extdata/organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID.gp.rds"))
###################
library(microtrait)
bioproject_id = ""
# isolates
#genomes_basedir = "/Users/ukaraoz/Work/microtrait/code/microtrait_extdata/genomic"
genomes_basedir = "/global/homes/u/ukaraoz/cscratch/alltarballs/organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID_taxonids"
genomes_dir = paste0(genomes_basedir, bioproject_id)
genomes_files = list.files(genomes_dir, full.names = T, recursive = F, pattern = "\\d+.fna$")
rds_dir = genomes_dir
genomes_files_done = file.path(genomes_dir, sub(".microtrait.rds", ".fna", list.files(rds_dir, pattern = "*rds$")))
genomes_files_todo = setdiff(genomes_files, genomes_files_done)
#genomes_files = genomes_files_remain
#results = run.parallel(genomes_files[1:2000], out_dir = genomes_dir, ncores = 10)
#indices = 8001:10306
run.parallel(genomes_files_todo, out_dir = dirname(genomes_files_todo), save_tempfiles = T, ncores = 35)

# MAGs
#/usr/bin/rm ./PRJNA362212_Dombrowski/genbank/*/*/*microtrait*
#/usr/bin/rm ./PRJNA362212_Dombrowski/genbank/*/*/*prodigal.faa
#/usr/bin/rm ./PRJNA362212_Dombrowski/genbank/*/*/*dbcan.domtblout
library(microtrait)
bioproject_id = "PRJNA348753_Parks"
# isolates
#genomes_basedir = "/Users/ukaraoz/Work/microtrait/code/microtrait_extdata/genomic"
genomes_basedir = "/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download"
genomes_dir = file.path(genomes_basedir, bioproject_id)
genomes_files = read.table(file.path(genomes_basedir, paste0(bioproject_id, ".fullpath_genomes.txt")))[,1]

genomes_files_done = sub(".microtrait.rds", ".fna", list.files(file.path(genomes_dir, "microtrait"), pattern = "*.rds$"))
genomes_files_todo = setdiff(basename(genomes_files), basename(genomes_files_done))
genomes_files_todo %in% basename(genomes_files)

genomes_files_remain = genomes_files[3500:length(genomes_files)]
run.parallel(fna_files = genomes_files_remain, out_dirs = rep(file.path(genomes_dir, "microtrait1"), length(genomes_files_remain)),
             save_tempfiles = T, ncores = 38)
extract.traits(genomes_files[1], out_dir = dirname(genomes_files)[1], save_tempfiles = T)

fna.dir = "/Users/ukaraoz/Work/microtrait/code/inst/extdata/genomic/in"
fna_files = list.files(fna.dir, full.names = T,pattern = "\\d+.fna$")
out_dirs = dirname(fna_files)[1:4]
run.parallel(fna_files = fna_files[1:4], out_dirs = dirname(fna_files)[1:4],
             save_tempfiles = T, ncores = 38)
parallel::mclapply(1:length(fna_files),
                   function(i) {
                     returnList = extract.traits(fna_files[i], out_dirs[i], save_tempfiles = save_tempfiles, type = type)
                     returnList
                   },
                   mc.cores = ncores)


run.parallel <- function(fna_files, out_dirs, save_tempfiles = F, type = "genomic", ncores) {
  tictoc::tic(paste0("Running microtrait for ", length(fna_files), " genomes"))
  parallel::mclapply(1:length(fna_files),
                     function(i) {
                       returnList = extract.traits(fna_files[i], out_dirs[i], save_tempfiles = save_tempfiles, type = type)
                       returnList
                     },
                     mc.cores = ncores)

  # alternative way to parallelize, more error-prone in cluster environment due to memory allocation errors
  #cl <- parallel::makeCluster(ncores)
  #doParallel::registerDoParallel(cl)
  #results = foreach(i = 1:length(fna_files), .packages = c("microtrait")) %dopar% {
  #  extract.traits(fna_files[i], save_tempfiles = T, out_dir = out_dir)
  #}
  #parallel::stopCluster(cl)
  tictoc::toc(log = "TRUE")
  #results
}

#rds_files = unlist(purrr::map_depth(results, 1, "rds_file"))
rds_files = list.files("/Users/ukaraoz/Work/microtrait/code/microtrait-analysis/img.rds", full.names = T, recursive = F, pattern = ".microtrait.rds$")
trait_matrix = combine.results(rds_files, type = "trait_atgranularity2")
hmm_matrix = combine.results(rds_files[1:10], type = "gene")

for(i in 1:length(genomes_files)) {
  r = extracttraits(genomes_files[i], type = "genomic", save_tempfiles = T, out_dir = cache_dir)
}
run.microtrait.parallel(genomes_files, out_dir = genomes_dir, ncores = 55)

library()
ncores = 4
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
results = foreach(i = 1:length(fna_files), .packages = c("microtrait")) %dopar% {
  extract.traits(fna_files[i], save_tempfiles = T, out_dir = cache_dir)
}

fna_files = genomes_files[1:4]
results = run.parallel(fna_files, out_dir = cache_dir, ncores = 4)


