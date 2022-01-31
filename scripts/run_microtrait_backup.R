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


