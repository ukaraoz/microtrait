# microTrait
microTrait is an R package that provides a workflow to extract fitness traits from microbial genome sequences. microTrait uses profile hidden Markov models (profile-HMM) and simple logical operations to predict and map protein family content represented in a genome sequence to fitness traits. The HMMs used in the microTrait framework ([microtrait-hmms](https://github.com/ukaraoz/microtrait-hmm)) represent protein family sequence diversity accumulated in genomes and metagenomes in [IMG/M](https://img.jgi.doe.gov/cgi-bin/m/main.cgi), have been curated, benchmarked using [KEGG](https://www.kegg.jp/) orthologs to establish *trusted cutoff (TC)* scores.

<img src="https://github.com/ukaraoz/microtrait/blob/master/microtrait.png" width="45%">

## Installation

### <a name="dependencies_overview"></a> Dependencies:
microTrait has the following software and data dependencies.
Software:

* **[HMMER](http://hmmer.org/) (>= v3.1b2)**
* **[Prodigal](https://github.com/hyattpd/Prodigal) (>= v2.6.3)**


These programs should be in your PATH so that they can be accessed regardless of location.

Data:

* **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)**: gene level profile-HMM database underlying microTrait framework
* **[dbCAN-HMMdb](http://bcb.unl.edu/dbCAN2/download/Databases/)**: domain level profile-HMM database for Carbohydrate-active enzymes.

See setup section below for deployment of these databases after installation.

### <a name="install"></a> Install:
The package will be available on CRAN soon. The developmental version can be installed using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package. To install the R package directly from GitHub, do in R:

```{r tidy = FALSE}
install.packages("devtools")
library(devtools)
devtools::install_github("ukaraoz/microtrait")
```

### <a name="setup"></a> Setup:
Due to its large size (~170M), **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)** isn't packaged with microTrait. After installation, **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)** database has to be downloaded and deployed. For this, we use **[piggyback](https://cran.r-project.org/web/packages/piggyback/index.html)** R package. Install and load **[piggyback](https://cran.r-project.org/web/packages/piggyback/index.html)**

```{r tidy = FALSE}
install.packages("piggyback")
library(piggyback)
```
microTrait includes a function (`prep.hmmmodels()`) that downloads and deploys hmm models that underlie microTrait. This function should be run once after the installation.

```{r tidy = FALSE}
prep.hmmmodels()
```
The following files should now be available under `microtrait/extdata/hmm` in your R install directory (`.libPaths()`):`dbCAN-HMMdb-V8.subset.hmm*, microtrait.hmm*`

## Usage
### Single genome
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

### Multiple (thousands) genomes (parallel workflow)

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
