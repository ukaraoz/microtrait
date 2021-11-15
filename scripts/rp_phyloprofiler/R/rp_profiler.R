library(microtrait)
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(parallel)

base = "/global/homes/u/ukaraoz/cscratch/alltarballs"
in.dir = file.path(base, "organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID_taxonids")
golddata = readRDS("/global/homes/u/ukaraoz/cscratch/alltarballs/organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID.metadata.rds")
# path to pullseq binary
pullseq = "/global/homes/u/ukaraoz/bin/pullseq/pullseq"

#source(file.path(base, "bin/rp_phyloprofiler/R/hmmer.R"))
source(file.path("/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/rp_phyloprofiler", "bin/R/utils.R"))
data.dir = "/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/rp_phyloprofiler/bin/data"

# path to hmm models
hmm.files = dir(file.path(data.dir, "hmm"), pattern = "hmm", full.names = T)
rp.loci = sub(".hmm", "", basename(hmm.files))
# path to the output directory
out.dir = file.path(base, paste0("organisms_habitatfiltered_matchedbyanalysis_wIMGTaxonID_taxonids", ".rp_phyloprofiler_out"))
if(!dir.exists(out.dir)) {
  dir.create(out.dir)
}

# path to the input files
faa.files = file.path(in.dir, dir(in.dir, recursive = F, pattern = ".genes.faa$"))
fna.files = file.path(in.dir, dir(in.dir, recursive = F, pattern = ".genes.fna$"))
#phylodist.files = file.path(base, dir(base, recursive = T, pattern = ".phylodist.txt$"))

# detect ribosomal protein loci
faa.files.detected = list()
faa.files.detected=parallel::mclapply(1:length(faa.files),
                     function(i) {
                        scan.faa(faa.files[i], hmm_files = grep("rpsC", hmm.files, value = T), out_dir = out.dir)
                     },
                     mc.cores = 20)
#faa.files.detected = lapply(faa.files, scan.faa, hmm_files = grep("rpsC", hmm.files, value = T), out_dir = out.dir)
faa.files.detected = do.call(bind_rows, faa.files.detected)
faa.files.detected1 = faa.files.detected %>% 
  dplyr::mutate(faa = stringr::str_replace(basename(faa), ".genes.faa", "")) %>%
  dplyr::arrange(faa) %>% dplyr::group_by(faa) %>% dplyr::slice(1) %>% #pick one copy per genome
  dplyr::left_join(golddata, by = c("faa" = "IMG Taxon ID")) %>%
  dplyr::mutate(header = paste0(`gene_name`, "_", `NCBI Superkingdom`, "_", `NCBI Phylum`, "_", 
                                `NCBI Class`, "_", `NCBI Order`, "_", `NCBI Family`, "_", 
                                `NCBI Genus`)) %>% 
  dplyr::ungroup() %>%
  dplyr::select(c(`gene_name`, `faa`, `header`)) %>%
  dplyr::mutate(pullseq = paste0("echo ", "\"", `gene_name`, "\"", 
    " |pullseq -i ", in.dir, "/", `faa`, ".genes.faa -N |sed 's/", `gene_name`, ".*/", `header`, "/'",
    ">", out.dir, "/", `faa`, ".rpsC.faa"))

pullseq_commands = faa.files.detected1 %>% dplyr::pull(pullseq)
write.table(pullseq_commands, 
                file.path(base, "pullseq_commands.rp_phyloprofiler_out.sh"),
                row.names = F, col.names = F, quote = F)

# MAGs
base = "/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download"
bioproject_ids = c("PRJNA288027_Anantharaman", "PRJNA348753_Parks", "PRJNA362212_Dombrowski", "PRJNA386568_Woodcroft")

for(i in 1:length(bioproject_ids)) {
  in.dir = file.path(base, bioproject_ids[i], "microtrait")
  out.dir = file.path(base, bioproject_ids[i], "rp_phyloprofiler_out")
  if(!dir.exists(out.dir)) {
    dir.create(out.dir)
  }
  faa.files = file.path(in.dir, dir(in.dir, recursive = F, pattern = ".prodigal.faa$"))
  # detect ribosomal protein loci
  faa.files.detected = list()
  faa.files.detected=parallel::mclapply(1:length(faa.files),
                       function(i) {
                          scan.faa(faa.files[i], hmm_files = grep("rpsC", hmm.files, value = T), out_dir = out.dir)
                       },
                       mc.cores = 20)
  faa.files.detected = do.call(bind_rows, faa.files.detected)
  faa.files.detected1 = faa.files.detected %>% 
    dplyr::mutate(faa = stringr::str_replace(basename(faa), "_genomic.prodigal.faa", "")) %>%
    dplyr::arrange(faa) %>% dplyr::group_by(faa) %>% dplyr::slice(1) %>% #pick one copy per genome
    dplyr::mutate(header = paste0(`gene_name`, "__", `faa`)) %>% 
    dplyr::ungroup() %>%
    dplyr::select(c(`gene_name`, `faa`, `header`)) %>%
    dplyr::mutate(pullseq = paste0("echo ", "\"", `gene_name`, "\"", 
      " |pullseq -i ", in.dir, "/", `faa`, "_genomic.prodigal.faa -N |sed 's/", `gene_name`, ".*/", `header`, "/'",
      ">", out.dir, "/", `faa`, ".rpsC.faa"))
  
  pullseq_commands = faa.files.detected1 %>% dplyr::pull(pullseq)
  write.table(pullseq_commands, 
                  file.path(base, bioproject_ids[i], "pullseq_commands.rp_phyloprofiler_out.sh"),
                  row.names = F, col.names = F, quote = F)
}





