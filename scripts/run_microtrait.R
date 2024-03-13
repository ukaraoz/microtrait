#!/usr/bin/env Rscript

# 
suppressPackageStartupMessages(library(argparser))
suppressWarnings(suppressMessages(library(microtrait)))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(parallel))
suppressMessages(library(vegan))
#suppressMessages(library(here))
#source(here::here("microtrait_kbaseutils.R"))
# create parser object

parser <- arg_parser(description = "Runs microtrait for genomes and defines guilds.", 
                     name = "run_microtrait",
                     hide.opts = TRUE)
parser <- add_argument(parser, "genome_fasta_files", type = "character", help = "input file of fasta paths for genomes")
parser <- add_argument(parser, "output_directory", type = "character", help = "output directory")
parser <- add_argument(parser, "dataset_name", type = "character", default = "microtrait", help = "dataset name")
parser <- add_argument(parser, "number_of_cores", type = "integer", default = 1, help = "number of cores")
#parser <- add_argument(parser, "--nguilds", type = "integer", help = "number of guilds")
#parser <- add_argument(parser, "--traitmatrix", type = "character", help = "precomputed trait matrix")
parser <- add_argument(parser, "--variance_explained", type = "integer", short = "ve", default = 70, help = "explained variance")
parser <- add_argument(parser, "--verbose", flag = T, short = "v", help = "be verbose")
args <- parse_args(parser)

genome_fasta_files = args$genome_fasta_files
output_directory = args$output_directory
dataset = args$dataset_name
ncores = args$number_of_cores
verbose = args$verbose
#traitmatrix_file = args$traitmatrix
variance_explained = args$variance_explained

genome_fasta_files = "/global/homes/u/ukaraoz/cscratch/Hopland39Isolates/Hopland39Isolates.fna.files.txt"
output_directory = "/global/homes/u/ukaraoz/cscratch/Hopland39Isolates/microtrait.out"
dataset = "rhizosphere_isolates"
ncores = 1
verbose = T
variance_explained = 70
#saveRDS(args, file = "~/Work/kbase/microtrait_wrapper/args.rds")
genome_files = read.table(genome_fasta_files)

genome_fasta_paths = scan(genome_fasta_files, sep="", what = "", quiet=TRUE)
genome_fasta_paths_ok = genome_fasta_paths[which(file.exists(genome_fasta_paths) == TRUE)]
if(verbose) {
  message(length(genome_fasta_paths_ok), " fasta files will be processed. ", 
          length(genome_fasta_paths)-length(genome_fasta_paths_ok), " didn't exist.")
}
microtrait_results = parallel::mclapply(1:length(genome_fasta_paths_ok),
                        function(i) {
                           returnList = extract.traits(genome_fasta_paths_ok[i],
                                                       output_directory,
                                                       save_tempfiles = F, 
                                                       type = "genomic", mode = "meta")
                           returnList
                        },
                        mc.cores = ncores)
rds_files = unlist(lapply(microtrait_results, "[[", "rds_file"))
genomeset_results = microtrait::make.genomeset.results(rds_files = rds_files,
                                                       ids = sub(".microtrait.rds", "", basename(rds_files)),
                                                       growthrate = T, optimumT = F,
                                                       ncores = ncores)
trait_matrixatgranularity3 = genomeset_results[["trait_matrixatgranularity3"]]
objects = c("trait_matrixatgranularity1", "trait_matrixatgranularity2", "trait_matrixatgranularity3", "hmm_matrix", "rule_matrix")
for(o in 1:length(objects)) {
  towrite = genomeset_results[[objects[o]]]
  write.table(towrite,
              file = file.path(output_directory, paste0(dataset, "", objects[o], ".txt")),
              row.names = F, col.names = T, sep = "\t", quote = F)
}
calc_intergenome_variance(genomeset_results = genomeset_results, outdir = output_directory)

