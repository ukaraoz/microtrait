# microTrait
microTrait is an R package that provides a workflow to extract fitness traits from microbial genome sequences. microTrait uses profile hidden Markov models (profile-HMM) and simple logical operations to predict and map protein family content represented in a genome sequence to fitness traits. The HMMs used in the microTrait framework ([microtrait-hmms](https://github.com/ukaraoz/microtrait-hmm)) represent protein family sequence diversity accumulated in genomes and metagenomes in [IMG/M](https://img.jgi.doe.gov/cgi-bin/m/main.cgi), have been curated, benchmarked using [KEGG](https://www.kegg.jp/) orthologs to establish *trusted cutoff (TC)* scores.

<img src="https://github.com/ukaraoz/microtrait/blob/master/microtrait.png" width="45%">

## Installation

### <a name="dependencies_overview"></a> Dependencies:
microTrait has the following software and data dependencies.

#### Software

The following binaries should be in your PATH so that they can be accessed regardless of location.

* **[HMMER](http://hmmer.org/) (>= v3.1b2)**
* **[Prodigal](https://github.com/hyattpd/Prodigal) (>= v2.6.3)**

The following are required for calculating features used in calculating Optimal Growth Temperatures based on regression models from [Sauer (2009)](https://doi.org/10.1093/bioinformatics/btz059)

* **[Infernal](http://eddylab.org/infernal/infernal-1.1.2.tar.gz) (= v1.1.2)**
* **[tRNAscan-SE](http://trna.ucsc.edu/software/trnascan-se-2.0.0.tar.gz) (= v2.0.0)**
* **[bedtools](https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz) (>= v2.27)**
* **[barrnap](https://github.com/tseemann/barrnap/archive/0.9.tar.gz) (= v0.9)**

#### Data
* **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)**: gene level profile-HMM database underlying microTrait framework
* **[dbCAN-HMMdb](http://bcb.unl.edu/dbCAN2/download/Databases/)**: domain level profile-HMM database for Carbohydrate-active enzymes.

See setup section below for deployment of these databases after installation.

### <a name="install"></a> Installation:
The developmental version can be installed using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package. To install the R package directly from GitHub, do in R:

```{r tidy = FALSE}
install.packages("devtools")
library(devtools)
devtools::install_github("ukaraoz/microtrait")
```

### <a name="setup"></a> Setup of HMM database:
Due to its large size (~170M), **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)** isn't packaged with microTrait. After installation, **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)** database has to be downloaded and deployed.

microTrait includes a function (`prep.hmmmodels()`) that downloads and deploys hmm models that underlie microTrait. This function should be run once after the installation.

```{r tidy = FALSE}
prep.hmmmodels()
```
The following HMM database files should now be available under `microtrait/extdata/hmm` in your default R library path (`.libPaths()`)

```{r tidy = FALSE}
list.files(file.path(.libPaths(), "microtrait/extdata/hmm/hmmpress"))
[1] "dbcan.select.v8.hmmdb.h3f" "dbcan.select.v8.hmmdb.h3i"
[3] "dbcan.select.v8.hmmdb.h3m" "dbcan.select.v8.hmmdb.h3p"
[5] "microtrait.hmmdb.h3f"      "microtrait.hmmdb.h3i"     
[7] "microtrait.hmmdb.h3m"      "microtrait.hmmdb.h3p"
```

## Usage
### Extracting traits from a genome sequence
#### Single genome
When microTrait is setup, you can extract traits from a genome sequence. For the following example, we use the genome sequence included as part of microTrait package.

```{r tidy = FALSE}
library(microtrait)
genome_file <- system.file("extdata/examples/2775507255/in/2775507255.fna", package="microtrait")
microtrait_result = extracttraits(genome_file)
```

microTrait packages its results into an R list with the following names components:

* `microtrait_result$id`: genome id, currently set to fasta filename for the genome sequence
* `microtrait_result$genes_detected`: microtrait-hmm models detected in the genome sequence (using the benchmarked TC thresholds)
* `microtrait_result$domains_detected`: CAZy domains detected in the genome sequence
* `microtrait_result$all_traits`: Extracted traits
* `microtrait_result$time_log`: Run time logs

#### Multiple (thousands) genomes (parallel workflow)

A simple data level parallelization of microtrait for multiple genomes can be easily implemented using `foreach` and `parallel` R packages. Using a parallel workflow, microTrait is able to process a few thousands of genomes within 2 hours on a 64 core Intel environment.

To exemplify the processing of multiple genomes with microtrait, we will use 2545 metagenome assembled genomes (MAGs) from [Anantharaman et al. 2016] (https://www.nature.com/articles/ncomms13219#ref20) ("Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system". The data has been deposited to NCBI under [PRJNA288027](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA288027). GenBank Assembly IDs are accessible in a downloadable table from [PRJNA288027](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA288027), and also provided [PRJNA288027-accession-list.txt](LINK).

To bulk download these genomes, we will use [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) as follows (replace genomes_dir with path to your output folder):

``
ncbi-genome-download -s genbank -F fasta,protein-fasta -A PRJNA386568-accession-list.txt --parallel 4 --retries 4 --output-folder genomes_dir all
gunzip -r $genomes_dir
``

First determine, on your machine/server, the number of "cores", or computational units with `detectCores()` function:

```{r detectCores, warning=FALSE, message=FALSE}
library("parallel")
detectCores()
```

```{r detectCores, warning=FALSE, message=FALSE}
library("tictoc")
#genomes_dir="/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download/test"
bioproject_id = "PRJNA288027"
genomes_dir = paste0("/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags/download/", bioproject_id)
genomes_files = list.files(genomes_dir, full.names = T, recursive = T, pattern = "genomic.fna$")
r=extracttraits(genomes_files[1])
cache_dir = "/global/homes/u/ukaraoz/m3260RR/microtrait-hmm/prep-ncbimags"
results <- mclapply(1:4,
             FUN=function(i) extracttraits(genomes_files[i]),
             mc.cores = 10)
message("Running on: ", system("cat /proc/cpuinfo | grep \"model name\" | uniq | sed 's/.*: //'", intern = TRUE), "\n")
message("Number of cores:", detectCores(), "\n")
tictoc::tic.clearlog()
tictoc::tic(paste0("Running microtrait for ", length(genomes_files), " genomes from ", bioproject_id))
result <- mclapply(1:length(genomes_files), function(i) {
        r = extracttraits(genomes_files[i])
        saveRDS(r, file = file.path(cache_dir, paste0(fs::path_file(genomes_files[i]), ".microtrait.rds")))
        r}, mc.cores = 62)
tictoc::toc(log = "TRUE")
Running microtrait for 2545 genomes from PRJNA288027: 4509.458 sec elapsed
}
```
### Building trait matrix from microtrait outputs
`microtrait::make.genomeset.results` combines microtrait outputs for multiple genomes into trait matrices (genomes x traits). Additionally, the resulting object will hold rules and hmm matrices. This specific example will use 80% of all cores available.

```{r make.genomeset.results, warning=FALSE, message=FALSE}
rds_files = list.files(file.path(base_dir, "microtrait-out"),
                       full.names = T, recursive = F, pattern = ".microtrait.rds$")
genomeset_results = make.genomeset.results(rds_files = rds_files,
                                           ids = sub(".microtrait.rds", "", basename(rds_files)),
                                           ncores = floor(0.8*detectCores()))
saveRDS(genomeset_results, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
```
If you have genome metadata, you can add them to the resulting matrices as additional columns using `microtrait::add.metadata`.

```{r add.metadata, warning=FALSE, message=FALSE}
genome_metadata = readRDS(file.path(base_dir, paste0(dataset, ".metadata.rds")))
genomeset_results_wmetadata = add.metadata(genomeset_results, genome_metadata, genome_metadata_idcol = "IMG Taxon ID") %>% convert_traitdatatype(binarytype = "logical")
saveRDS(genomeset_results_wmetadata, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
```
### Analyzing trait data from a set of genomes
#### Normalize count traits
The utility function `microtrait::trait.normalize` will correct the count traits for a given metric from the metadata column. Here, we use genome size (*Estimated Size*) to normalize count traits.

```{r trait.normalize, warning=FALSE, message=FALSE}
genomeset_results_wmetadata_norm = genomeset_results_wmetadata %>% trait.normalize(normby = "Estimated Size")
```
#### Trait-to-trait associations
Given the trait values across a large number of genomes, we can explore the association structure between traits. For this analysis, we use the trait measurements at the most granular level (*trait_matrixatgranularity3*) and focus on non-aquatic genomes.

```{r trait2traitcorrelations, warning=FALSE, message=FALSE}
trait_matrixatgranularity3 = genomeset_results_wmetadata_norm[["trait_matrixatgranularity3"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  dplyr::select(-`Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation`) %>%
  dplyr::select(c(1:189))
```

We first transform continous traits into binary traits, treating them as presence/absence traits.

```{r trait2traitcorrelations, warning=FALSE, message=FALSE}
trait_matrixatgranularity3_binary = trait_matrixatgranularity3 %>% trait.continuous2binary
```

Use `microtrait::trait2trait_corr` to plot trait correlation matrix and write correlations into a tab-delimited file.

```{r trait2traitcorrelations, warning=FALSE, message=FALSE}
microtrait::trait2trait_corr(trait_matrixatgranularity3_binary, verbose = TRUE, idcol = "id", outdir = base_dir, dataset = "soilgenomes")
```

#### Defining guilds
##### Generate distance matrix
Use `microtrait::calc_mixeddist` to generate the genome to genome distance matrix based on the trait profiles. for The default distance metric is [Wishart (2003)](https://doi.org/10.1007/978-3-642-55721-7_23) distance.

```
genomeset_distances = trait_matrixatgranularity3 %>% microtrait::calc_mixeddist()
```

```
clusters_traitmatrix = cluster_traitmatrix(trait_matrixatgranularity3_binary,
                                  clustering_distance_rows = genomeset_distances,
                                  clustering_distance_cols = "binary",
                                  outdir = base_dir, dataset = "soilgenomes")
```

We quantify the variance in the distance matrix as a function of the number of guilds. This is compute-heavy, a compute environment with high parallelization capability is warrented.
First generate a range for number of guilds across which we want to quantify variance. As an example, we go from 2 guilds to the total number of genomes with a step size of 2.

```
nguilds = seq(2, nrow(trait_matrixatgranularity3), 2)
hclust_rows = clusters_traitmatrix[[1]]
maov_results_list = parallel::mclapply(2:length(nguilds),
                     function(i) {
                       v = cutree(hclust_rows, nguilds[i])
                       genome2guild = data.frame(guild = factor(v))
                       rownames(genome2guild) = names(v)
                       adonis_results = vegan::adonis2(distance ~ guild, data = genome2guild, perm = 1)
                       adonis_results
                     },
                     mc.cores = floor(0.8*detectCores()))
maov_results = matrix(nrow = 0, ncol = 4)
for(i in 1:length(maov_results_list)) {
  #cat(i, "\n")
  maov_results = rbind(maov_results,
                       c(nguilds[i], maov_results_list[[i]]$R2))
}

maov_results = data.frame(`number of guilds` = maov_results[,1],
                          `percent variance` = round(maov_results[,2]*100,2),
                          check.names = F)
p = ggplot2::ggplot(maov_results,aes(x=`number of guilds`, y=`percent variance explained`)) +
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
dataset = "soilgenomes"
ggsave(p, height = 6, width = 12,
       filename = file.path(base_dir, paste0(dataset, ".guilds.percentvariance.pdf")))
```

Use `microtrait::defineguilds` to define guilds at a given percent variance explained. We will define guilds so as to capture 70% of inter-genome trait variation.

```
nguildsat70perc = maov_results[which.min(abs(maov_results[,2]-70)), "number of guilds"]
defined_guilds = microtrait::define_guilds(trait_matrixatgranularity3_binary,
                                        clusters_traitmatrix,
                                        nguilds = nguildsat70perc,
                                        guildsizecutoff = 50,
                                        outdir, dataset, verbose = TRUE)
```
