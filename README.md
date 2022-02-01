# *microTrait*
*microTrait* is an R package that provides a workflow to extract fitness traits from microbial genome sequences. *microTrait* uses profile hidden Markov models (profile-HMM) and simple logical operations to predict and map protein family content represented in a genome sequence to fitness traits. The HMMs used in the *microTrait* framework ([microtrait-hmms](https://github.com/ukaraoz/microtrait-hmm)) represent protein family sequence diversity accumulated in genomes and metagenomes in [IMG/M](https://img.jgi.doe.gov/cgi-bin/m/main.cgi), have been curated, benchmarked using [KEGG](https://www.kegg.jp/) orthologs to establish *trusted cutoff (TC)* scores.

The genome inferred traits represented by *microTrait* are summarized in the figure below. 

<img src="https://github.com/ukaraoz/microtrait/blob/master/microtrait.png" width="60%">

*microTrait* maps a genome sequence to a hierarchy of traits summarized in this [figure](https://github.com/ukaraoz/microtrait/blob/master/microtrait-hierarchy.pdf). The genomic basis for each trait was compiled from the literature accessible in [this file](https://github.com/ukaraoz/microtrait/blob/master/inst/extdata/ST8.microtrait_references.txt).



## Installation

*microTrait* has the following software and data dependencies.

### 1. Software dependencies

The following binaries should be in your PATH so that they can be accessed regardless of location.

* **[HMMER](http://hmmer.org/) (>= v3.1b2)**
* **[Prodigal](https://github.com/hyattpd/Prodigal) (>= v2.6.3)**

The following are required for calculating features used in calculating Optimal Growth Temperatures based on regression models from [Sauer (2009)](https://doi.org/10.1093/bioinformatics/btz059)

* **[Infernal](http://eddylab.org/infernal/infernal-1.1.2.tar.gz) (= v1.1.2)**
* **[tRNAscan-SE](http://trna.ucsc.edu/software/trnascan-se-2.0.0.tar.gz) (= v2.0.0)**
* **[bedtools](https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz) (>= v2.27)**
* **[barrnap](https://github.com/tseemann/barrnap/archive/0.9.tar.gz) (= v0.9)**

### 2. Data dependencies

* **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)**: gene level profile-HMM database underlying microTrait framework
* **[dbCAN-HMMdb](http://bcb.unl.edu/dbCAN2/download/Databases/)**: domain level profile-HMM database for Carbohydrate-active enzymes.

See setup section below for deployment of these databases after installation.

### 3. R package dependencies

The developmental version of *microtrait* can be installed using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package. First install and load [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

```{r tidy = FALSE}
install.packages("devtools")
library(devtools)
```
Next install R package dependencies for *microTrait*:

*  [CRAN](https://cran.r-project.org/) package dependencies

	Missing dependencies available on [CRAN](https://cran.r-project.org/) can be installed as follows:

	```{r tidy=TRUE}
	list_of_packages = c("R.utils", "RColorBrewer", "ape", "assertthat", "checkmate",
	 "coRdon", "corrplot", "doParallel", "dplyr", "futile.logger", "grid", "gtools", "kmed", "lazyeval", "magrittr", "parallel", "pheatmap", "readr", "stringr", "tibble", 
	  "tictoc", "tidyr")
	newpackages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) install.packages(newpackages)```

* [Bioconductor](https://www.bioconductor.org/) package dependencies

	*microTrait* also depends on [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [coRdon](https://github.com/BioinfoHR/coRdon), and [ComplexHeatmap](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html), which are [Bioconductor](https://bioconductor.org/) packages. To install them, run the following:

	```{r tidy = FALSE}
	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")
	BiocManager::install("Biostrings")
	BiocManager::install("coRdon")
	BiocManager::install("ComplexHeatmap")
	```

* Developmental package dependencies

	*microTrait* uses [gRodon](https://github.com/jlw-ecoevo/gRodon) for estimation of maximum growth rates from genomic sequences. Install with [devtools](https://github.com/r-lib/devtools):

	```{r tidy = FALSE}
	devtools::install_github("jlw-ecoevo/gRodon")
	```

### 4. Installation of *microTrait*:

Now install *microTrait* with [devtools](https://github.com/r-lib/devtools):

```{r tidy = FALSE}
devtools::install_github("ukaraoz/microtrait")
```

### 5. <a name="setup"></a> Setup of HMM database:
Due to its large size (~170M), **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)** isn't packaged with *microTrait*. After installation, **[microtrait-hmm](https://github.com/ukaraoz/microtrait-hmm)** database has to be downloaded and deployed.

Load *microTrait*:

```{r tidy = FALSE}
library(microtrait)
```

*microTrait* includes a function (`prep.hmmmodels()`) that downloads and deploys hmm models that underlie *microTrait*. This function should be run once after the installation.

```{r tidy = FALSE}
microtrait::prep.hmmmodels()
```
The following HMM database files should now be available under `microtrait/extdata/hmm` in your default R library path (`.libPaths()`)

```{r tidy = TRUE}
list.files(file.path(.libPaths(), "microtrait/extdata/hmm/hmmpress"))
 [1] "arcbacribosomal.hmmdb.h3f" "arcbacribosomal.hmmdb.h3i"
 [3] "arcbacribosomal.hmmdb.h3m" "arcbacribosomal.hmmdb.h3p"
 [5] "dbcan.select.v8.hmmdb.h3f" "dbcan.select.v8.hmmdb.h3i"
 [7] "dbcan.select.v8.hmmdb.h3m" "dbcan.select.v8.hmmdb.h3p"
 [9] "microtrait.hmmdb.h3f"      "microtrait.hmmdb.h3i"     
[11] "microtrait.hmmdb.h3m"      "microtrait.hmmdb.h3p"     
```

## Usage
### Extracting traits from a genome sequence
#### Single genome
Once *microTrait* is correctly setup, you can extract traits from a genome sequence. For the following example, we use the genome sequence of "Achromobacter xylosoxidans A-8 (IMG Taxon OID 2695420375)" included as part of *microTrait* package.

```{r tidy = FALSE}
library(microtrait)
genome_file <- system.file("extdata/genomic/2695420375.fna", package="microtrait")
microtrait_result = extract.traits(genome_file)
```

*microTrait* packages its results into an R list with the following names components:

* `microtrait_result$id`: genome id, currently set to fasta filename for the genome sequence
* `microtrait_result$genes_detected`: microtrait-hmm models detected in the genome sequence (using the benchmarked TC thresholds)

* `microtrait_result$genes_detected`: microtrait-hmm models detected in the genome sequence (using the benchmarked TC thresholds)
* `microtrait_result$domains_detected`: CAZy domains detected in the genome sequence
* `microtrait_result$all_traits`: Extracted traits
* `microtrait_result$time_log`: Run time logs

#### Multiple (thousands) genomes (parallel workflow)

A simple data level parallelization of *microTrait* for multiple genomes is implemented using `foreach` and `parallel` R packages in the `microtrait::extract.traits.parallel`` function.

<!--Using a parallel workflow, *microTrait* is able to process a few thousands of genomes within 2 hours on a 64 core Intel environment.-->


To exemplify the processing of multiple genomes with *microTrait*, we will process a test dataset of 100 genomes from IMG available [here](https://100genomes.s3.us-west-1.amazonaws.com/100genomes.tar.gz). Download and unpack it under your package install directory:


```{r download100genomes, warning=FALSE, message=FALSE, tidy = TRUE}
download.file("https://github.com/ukaraoz/microtrait-hmm/releases/download/latest/100genomes.tar.gz",
    system.file("extdata/genomic/100genomes.tar.gz", package = "microtrait"))
untar(system.file("extdata/genomic/100genomes.tar.gz", package = "microtrait"), 
      exdir = system.file("extdata/genomic", package = "microtrait"))
```

The downloaded files will be under the following `genomes_dir` directory:

```{r}
genomes_dir = system.file("extdata/genomic/100genomes", package = "microtrait")
genomes_files = list.files(genomes_dir, full.names = T, recursive = T, pattern = ".fna$")
```
	
<!--
To exemplify the processing of multiple genomes with microtrait, we will use 2545 metagenome assembled genomes (MAGs) from [Anantharaman et al. 2016] (https://www.nature.com/articles/ncomms13219#ref20) ("Thousands of microbial genomes shed light on interconnected biogeochemical processes in an aquifer system". The data has been deposited to NCBI under [PRJNA288027](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA288027). GenBank Assembly IDs are accessible in a downloadable table from [PRJNA288027](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA288027), and also provided [PRJNA288027-accession-list.txt](LINK).

To bulk download these genomes, we will use [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) as follows (replace genomes_dir with path to your output folder):

``
ncbi-genome-download -s genbank -F fasta,protein-fasta -A PRJNA386568-accession-list.txt --parallel 4 --retries 4 --output-folder genomes_dir all
gunzip -r $genomes_dir
-->

Determine the number of "cores", or computational units on your machine/server using [`parallel::detectCores()`](https://www.rdocumentation.org/packages/parallel/versions/3.6.2/topics/detectCores) function:

```{r detectCores, warning=FALSE, message=FALSE}
message("Number of cores:", parallel::detectCores(), "\n")
```

Here, we use ~70% of the available cores to process 100 genomes:

```{r extract.traits.parallel, warning=FALSE, message=FALSE}
library("tictoc")
tictoc::tic.clearlog()
tictoc::tic(paste0("Running microtrait for ", length(genomes_files)))
microtrait_results = extract.traits.parallel(genomes_files, dirname(genomes_files), ncores = floor(parallel::detectCores()*0.7)
tictoc::toc(log = "TRUE")
```

When completed, `microtrait_results` will hold a list of lists for *microTrait* results. We can pull the paths to the corresponding "rds" files as follows:

<!--
rds_files = list.files(system.file("extdata/genomic/100genomes",package = "microtrait"),
                       full.names = T, recursive = F, pattern = ".microtrait.rds$")
-->


```{r, tidy = TRUE}
rds_files = unlist(mclapply(microtrait_results, "[[", "rds_file", mc.cores = ncores))
```

Next, using these `rds_files` we will combine these results from multiple genomes.

### Building trait matrices from *microTrait* outputs
`microtrait::make.genomeset.results()` combines microtrait outputs for multiple genomes into trait matrices (genomes x traits). Additionally, the resulting object will hold rules and hmm matrices.

```{r make.genomeset.results, warning=FALSE, message=FALSE}
genomeset_results = make.genomeset.results(rds_files = rds_files,
                                           ids = sub(".microtrait.rds", "", basename(rds_files)),
                                           ncores = 1)
```

The output `genomeset_results` will be a list with 5 elements, each of which is data.frame corresponding to extracted traits at different granularities, rules, and microtrait-hmm hits.

```{r warning=FALSE, message=FALSE}
names(genomeset_results)
[1] "trait_matrixatgranularity1" "trait_matrixatgranularity2"
[3] "trait_matrixatgranularity3" "hmm_matrix"                
[5] "rule_matrix"
```

------
### Analyzing trait data from a set of genomes
To demonstrate analysis of trait data from a set of genomes, we will use a dataset consisting of all environmental isolate genomes from IMG database. The process to combine *microTrait* results for multiple genomes is as above. Here, we use the precomputed `microTrait` outputs for this set that are available as part of the `microTrait` package.

```{r make.genomeset.results, warning=FALSE, message=FALSE}
base_dir = system.file("extdata/genomic/precomputed",package = "microtrait")
dataset = "environmentalgenomes"
genomeset_results = readRDS(file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
lapply(genomeset_results, dim)
```

If you have genome metadata, you can append them to the resulting matrices as additional columns using `microtrait::add.metadata`. Here we use GOLD metadata for these genomes, also available as part of the package.

```{r add.metadata, warning=FALSE, message=FALSE}
library(dplyr)
genome_metadata = readRDS(file.path(base_dir, paste0(dataset, ".metadata.rds")))
genomeset_results_wmetadata = microtrait::add.metadata(genomeset_results, genome_metadata, genome_metadata_idcol = "IMG Taxon ID") %>% convert_traitdatatype(binarytype = "logical")
saveRDS(genomeset_results_wmetadata, file.path(base_dir, paste0(dataset, ".microtraitresults.rds")))
```
#### Normalize count traits
The utility function `microtrait::trait.normalize()` will correct the count traits for a given metric from the metadata column. Here, we use genome size (*Estimated Size*) to normalize count traits.

```{r trait.normalize, warning=FALSE, message=FALSE}
genomeset_results_wmetadata_norm = genomeset_results_wmetadata %>% trait.normalize(normby = "Estimated Size")
```
#### Trait-to-trait associations
Given the trait values across a large number of genomes, we can explore the association structure between traits. 

For this analysis, we use the trait measurements at the most granular level (*trait_matrixatgranularity3*) and focus on soil genomes.

```{r trait2traitcorrelations, warning=FALSE, message=FALSE}
traits = traits_listbygranularity[[3]] %>%
  dplyr::select(`microtrait_trait-name`) %>%
  dplyr::filter(`microtrait_trait-name` != "Resource Use:Chemotrophy:chemolithoautotrophy:anaerobic ammonia oxidation") %>%
  dplyr::pull(`microtrait_trait-name`)
trait_matrixatgranularity3 = genomeset_results_wmetadata_norm[["trait_matrixatgranularity3"]] %>%
  dplyr::filter(`Ecosystem_Category` %in% c("Terrestrial", "Plants")) %>%
  dplyr::filter(`Ecosystem_Type` %in% c("Soil", "Rhizoplane", "Rhizosphere", "Roots")) %>%
  dplyr::select(c("id", traits, "mingentime", "ogt")) %>%
  #dplyr::slice(1:100) %>%
  dplyr::filter(`mingentime` < 100) %>%   # max out mingentime at 25 days to avoid errors
  dplyr::filter(!is.na(`mingentime`) & !is.na(`ogt`))
```

We first transform continous traits into binary traits (with `microtrait:: trait.continuous2binary()`), treating them as presence/absence traits.

```{r trait2traitcorrelations, warning=FALSE, message=FALSE}
trait_matrixatgranularity3_binary = trait_matrixatgranularity3 %>% microtrait::trait.continuous2binary()

```

Use `microtrait::trait2trait_corr()` to plot trait correlation matrix and write correlations into a tab-delimited file.

```{r trait2traitcorrelations, warning=FALSE, message=FALSE}
trait2trait_corr(trait_matrixatgranularity3_binary, verbose = TRUE, idcol = "id", outdir = base_dir, dataset = "soilgenomes")
```

<img src="https://github.com/ukaraoz/microtrait/blob/master/trait2trait-correlation.png" width="60%">

#### Defining guilds
##### Generate genome-to-genome distance matrix
Use `microtrait::calc_mixeddist()` to generate the genome to genome distance matrix based on the trait profiles. The default distance metric is [Wishart (2003)](https://doi.org/10.1007/978-3-642-55721-7_23) distance for mixed data types.

```
genomeset_distances = trait_matrixatgranularity3 %>% microtrait::calc_mixeddist(idcol = "id", col2ignore = c("mingentime", "ogt"), method = "wishart", binarytype = "logical", byrow = 1, verbose = TRUE)
```

Also compute trait prevalence in the dataset with `microtrait::compute.prevalence()`:
```
prevalence = compute.prevalence(trait_matrixatgranularity3_binary, type="trait_matrixatgranularity3")
```

Now, apply hierarchical clustering based on the distance matrix and plot the corresponding heatmap. We use `unit` function from `grid` package.

```
library(grid)
library(ComplexHeatmap)
library(tibble)
A4_ratio = 11.75/8.25
width = grid::unit(8.25*4, "inches")
height = grid::unit(8.25*4*A4_ratio, "inches")
cluster_traitmatrix(trait_matrix = trait_matrixatgranularity3_binary,
                    idcol = "id", annot_cols = c("mingentime", "ogt"), granularity = 3,
                    clustering_distance_rows = genomeset_distances, clustering_distance_cols = "binary",
                    width = width, height = height,
                    heatmap_width = width*0.8, heatmap_height = height*0.70,
                    row_dend_width = width*0.05, column_dend_height = height*0.02,
                    rightannotation_width = width*0.1, topannotation_height = height*0.02,
                    bottomannotation_height = height*0.005,
                    outdir = base_dir, dataset = "soilgenomes", pdf = TRUE)
```
<img src="https://github.com/ukaraoz/microtrait/blob/master/traitmatrix.png" width="100%">


##### Quantification of inter-guild variance

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
```

We will read precomputed results:

```
maov_results_list = readRDS(system.file("extdata/genomic/precomputed/soilgenomes_maov_results_list.rds", package="microtrait"))

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
<img src="https://github.com/ukaraoz/microtrait/blob/master/traitvariance.png" width="60%">

Use `microtrait::defineguilds()` to define guilds at a given percent variance explained. We will define guilds so as to capture 70% of inter-genome trait variation.

```
nguildsat70perc = maov_results[which.min(abs(maov_results[,2]-70)), "number of guilds"]
defined_guilds = microtrait::define_guilds(trait_matrixatgranularity3_binary,
                                        clusters_traitmatrix,
                                        nguilds = nguildsat70perc,
                                        guildsizecutoff = 50,
                                        outdir, dataset, verbose = TRUE)
```

<img src="https://github.com/ukaraoz/microtrait/blob/master/guildsize-dist.png" width="100%">

<img src="https://github.com/ukaraoz/microtrait/blob/master/guildprofiles.png" width="100%">

