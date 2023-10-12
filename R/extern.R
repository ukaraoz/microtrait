#' @name run.hmmsearch
#'
#' @title Run hmmsearch (HMMER 3)
#'
#' @description Takes a fasta file and a Hidden Markov Model profile and
#' performs a search of the former over the latter.
#'
#' @param fasta A protein fasta file.
#'
#' @param hmm A hmm file. Must be pressed (see hmmpress from HMMER manual).
#'
#' @param domtblout_file
#'
#' @param n_threads An \code{integer}. The number of cores to use.
#' @importFrom assertthat assert_that
#' @return The path to a temporary file where the hmmsearch output is placed.
#'
#' hmmsearch --cpu $ncpus --domtblout $faafile".microtrait.domtblout" $hmm $faafile > /dev/null
run.hmmsearch <- function(faa_file = system.file("extdata", "2619619645.genes.faa", package = "microtrait", mustWork = TRUE),
                          hmm = NULL,
                          domtblout_file = NULL,
                          n_threads = 1) {
  # Check arguments
  faa_file = tryCatch({
    assertthat::assert_that(file.exists(faa_file))
    faa_file
  },
  error = function(e) {
    message("Protein fasta file ", `faa_file`, " doesn't exist.")
    print(e)
  }
  )

  hmm = match.arg(hmm, c("microtrait", "dbcan", "ribosomal"))
  hmm_file = switch(hmm,
                    #"microtrait" = file.path("/Users/ukaraoz/Work/microtrait/code/microtrait", "inst/extdata/hmm/hmmpress/microtrait.hmmdb"),
                    #"dbcan" = file.path("/Users/ukaraoz/Work/microtrait/code/microtrait", "inst/extdata/hmm/hmmpress/dbcan.select.v8.hmmdb")
                    "microtrait" = file.path(system.file(package="microtrait"), "extdata/hmm/hmmpress/microtrait.hmmdb"),
                    "dbcan" = file.path(system.file(package="microtrait"), "extdata/hmm/hmmpress/dbcan.select.v8.hmmdb"),
                    "ribosomal" = file.path(system.file(package="microtrait"), "extdata/hmm/hmmpress/arcbacribosomal.hmmdb")
  )
  #hmmsearch --cpu $ncpus --domtblout $faafile".microtrait.domtblout" $hmm $faafile > /dev/null
  result <- NULL
  if(available.external("hmmsearch")) {
    if(is.null(domtblout_file)) domtblout_file <- tempfile(pattern = "microtrait", fileext = ".domtblout")
    #command <- paste("hmmsearch",
    #                 "--cpu", n_threads,
    #                 "--domtblout", domtblout_file,
    #                 hmm,
    #                 faa_file,
    #                 ">/dev/null")
    message("Running hmmsearch for ", fs::path_file(faa_file), " with ", hmm, " models.")
    #message(command)

    #if(hmm == "microtrait") {
    #  system2("hmmsearch",
    #          args = c("--domtblout", domtblout_file,
    #                   "--cpu", n_threads,
    #                   "--noali",
    #                   "--cut_tc ", # CRUCIAL, otherwise doesn't use model specific thresholds
    #                   hmm_file, faa_file
    #          ), stdout = "/dev/null")
    #} else if(hmm == "dbcan") {
    #  system2("hmmsearch",
    #          args = c("--domtblout", domtblout_file,
    #                   "--cpu", n_threads,
    #                   "--noali",
    #                   hmm_file, faa_file
    #          ), stdout = "/dev/null")
    #}
    #file.remove(tmp.file)
    switch(hmm,
           microtrait = system2("hmmsearch",
                                args = c("--domtblout", domtblout_file,
                                         "--cpu", n_threads,
                                         "--noali",
                                         "--cut_tc ", # CRUCIAL, otherwise doesn't use model specific thresholds
                                         hmm_file, faa_file
                                ), stdout = "/dev/null"),
           dbcan = system2("hmmsearch",
                           args = c("--domtblout", domtblout_file,
                                    "--cpu", n_threads,
                                    "--noali",
                                    hmm_file, faa_file
                           ), stdout = "/dev/null"),
           ribosomal = system2("hmmsearch",
                               args = c("--domtblout", domtblout_file,
                                        "--cpu", n_threads,
                                        "--noali",
                                        "--cut_ga ",
                                        hmm_file, faa_file
                               ), stdout = "/dev/null"),
           stop("Invalid `hmm` value")
    )
    return(domtblout_file)
  }
}

#' @name run.prodigal
#' @title Finding coding genes
#'
#' @description Finding coding genes in genomic DNA using the Prodigal software.
#'
#' @param genome.file A FASTA file with the genome sequence(s).
#' @param faa.file If provided, prodigal will output all proteins to this fasta-file (text).
#' @param proc Either \code{"single"} or \code{"meta"}, see below.
#' @param mask.N Turn on masking of N's (logical)
#' @param bypass.SD Bypass Shine-Dalgarno filter (logical)
#'
#' @details The external software Prodigal is used to scan through a prokaryotic genome to detect the protein
#' coding genes. This free software can be installed from https://github.com/hyattpd/Prodigal.
#'
#' In addition to the standard output from this function, FASTA files with protein and/or DNA sequences may
#' be produced directly by providing filenames in \code{faa.file} and \code{ffn.file}.
#'
#' The input \code{proc} allows you to specify if the input data should be treated as a single genome
#' (default) or as a metagenome.
#'
#' The translation table is by default 11 (the standard code), but table 4 should be used for Mycoplasma etc.
#'
#' The \code{mask.N} will prevent genes having runs of N inside. The \code{bypass.SD} turn off the search
#' for a Shine-Dalgarno motif.
#'
#' @importFrom assertthat assert_that
#' @returns prodigal outputs
#'
#' @note The prodigal software must be installed on the system for this function to work, i.e. the command
#' \samp{system("prodigal -h")} must be recognized as a valid command if you run it in the Console window.
#'
run.prodigal <- function(genome_file = system.file("extdata/examples/2619619645/in", "2619619645.genes.fna", package = "microtrait", mustWork = TRUE),
                         fa_file = gsub(".fna", ".prodigal.fa", genome_file),
                         faa_file = gsub(".fna", ".prodigal.faa", genome_file),
                         mode = "single",
                         transtab = 11,
                         maskN = FALSE,
                         bypassSD = FALSE) {
  # Check arguments
  genome_file = tryCatch({
    assertthat::assert_that(file.exists(genome_file))
    genome_file
  },
  error = function(e) {
    message("Genome file ", `genome_file`, " doesn't exist.")
    print(e)
  }
  )

  faa_file = tryCatch({
    assertthat::assert_that(fs::dir_exists(fs::path_dir(faa_file)))
    faa_file
  },
  error = function(e) {
    message("Output file path ", fs::path_dir(faa_file), " isn't a valid path.")
    print(e)
  }
  )

  mode = tryCatch({
    assertthat::assert_that(mode %in% c("single", "meta"))
    mode
  },
  error = function(e) {
    message("Mode has to be either `single` or `meta`.")
    print(e)
  }
  )
  #checkmate::qexpect(genome_file, "S1", info = "info", label = "Genome fasta file")
  #checkmate::qassert(faa_file, "S1")
  #checkmate::assertChoice(mode, c("single", "meta"))
  #checkmate::qassert(maskN, "B")
  #checkmate::qassert(bypassSD, "B")
  #n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", parallel::detectCores(), "]")))
  #prodigal -i ./2775507255/2775507255.fna -o ./2775507255.prodigal -a ./2775507255.prodigal.faa -p single -q -f gff -m -n -g 11

  if(available.external("prodigal")){
    maskN <- ifelse(maskN, "-m", "")
    bypassSD <- ifelse(bypassSD, "-n", "")
    message("Running prodigal for ", fs::path_file(genome_file))
    # to do: reimplement with processx
    system2("prodigal",
            args = c("-f", "gff",
                     "-q",
                     "-g", transtab,
                     "-p", mode,
                     maskN, bypassSD,
                     "-i", genome_file,
                     "-d", fa_file,
                     "-a", faa_file,
                     "-o", tempfile(pattern = "prodigal", fileext = "gff")
            ),
            stdout = "/dev/null",stderr = "/dev/null"
    )
    #file.remove(tmp.file)
    return(list(fa_file = fa_file, faa_file = faa_file))
  }
}


run.tRNAscan <-function(genome_file) {
  # Check arguments
  genome_file = tryCatch({
    assertthat::assert_that(file.exists(genome_file))
    genome_file
  },
  error = function(e) {
    message("Genome file ", `genome_file`, " doesn't exist.")
    print(e)
  }
  )

  result <- NULL
  if(available.external("tRNAscan") & available.external("cmsearch")){
    message("Running tRNAscan-SE for ", fs::path_file(genome_file))
    # to do: reimplement with processx
    tempfilebase = tempfile(pattern = "tRNAscan")
    tempfile_trnascanstr = paste0(tempfilebase, ".str")
    tempfile_trnascangff = paste0(tempfilebase, ".gff")
    tempfile_trnascanfa = paste0(tempfilebase, ".fa")
    # final output with lowercase nucleotides kept as they are
    # pulled and converted to fasta from .str
    tempfile_trnascanlcfa = paste0(tempfilebase, ".lc.fa")
    tempfile_trnascanstderr = paste0(tempfilebase, ".stderr")
    #genome_file = "/opt/OGT_prediction-1.0.2/prediction/genomes/cyanobacterium_stanieri/2503283023.fa"
    system2("tRNAscan-SE",
            args = c("-G --brief -q ",
                     "--struct", tempfile_trnascanstr,
                     "--output", tempfile_trnascangff,
                     "--fasta", tempfile_trnascanfa,
                     " ", genome_file),
            stdin = "",
            stdout = "/dev/null",
            stderr = tempfile_trnascanstderr
    )
    check = grep("No tRNAs found", readLines(tempfile_trnascanstderr))
    if(length(check) == 0) {
      found = as.logical(TRUE)
    } else {
      found = as.logical(FALSE)
    }
    temp = readLines(tempfile_trnascanstr)
    temp = temp[grep('^Seq:.*', temp)]
    towrite = c(rbind(paste0(">tRNA_", 1:length(temp)), sub("Seq: ", "", temp)))
    write(x = towrite, file = tempfile_trnascanlcfa)

    #file.remove(tmp.file)
    result = list(found = found,
                  tRNA_str = tempfile_trnascanstr,
                  tRNA_gff = tempfile_trnascangff,
                  tRNA_fa = tempfile_trnascanfa,
                  tRNA_lcfa = tempfile_trnascanlcfa)
    return(result)
  }
}


#barrnap --quiet --threads 1 --kingdom bac /Users/ukaraoz/Work/microtrait/code/github/test/microtrait/inst/extdata/genomic/2503283023.fna \
#> ./output/genomes/desulfotomaculum_kuznetsovii/Desulfotomaculum_kuznetsovii_dsm_6115_ASM21470v1_dna_toplevel/barrnap/Bacteria.txt
run.barrnap <- function(genome_file = system.file("extdata", "2503283023.fna", package = "microtrait", mustWork = TRUE),
                        threads = 1,
                        kingdom = "bac",
                        quiet = ON) {
  # Check arguments
  genome_file = tryCatch({
    assertthat::assert_that(file.exists(genome_file))
    genome_file
  },
  error = function(e) {
    message("Genome file ", `genome_file`, " doesn't exist.")
    print(e)
  }
  )

  kingdom = tryCatch({
    assertthat::assert_that(kingdom %in% c("euk", "bac", "arc", "mito"))
    kingdom
  },
  error = function(e) {
    message("Kingdom has to be one of `euk`, `bac`, arc`, `mito` .")
    print(e)
  }
  )
  #checkmate::qexpect(genome_file, "S1", info = "info", label = "Genome fasta file")
  #checkmate::qassert(faa_file, "S1")
  #checkmate::assertChoice(mode, c("single", "meta"))
  #checkmate::qassert(maskN, "B")
  #checkmate::qassert(bypassSD, "B")
  #n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", parallel::detectCores(), "]")))
  #prodigal -i ./2775507255/2775507255.fna -o ./2775507255.prodigal -a ./2775507255.prodigal.faa -p single -q -f gff -m -n -g 11

  result <- NULL
  if(available.external("barrnap") & available.external("bedtools")){
    message("Running barrnap for ", fs::path_file(genome_file))
    # to do: reimplement with processx
    tempfilebase = tempfile(pattern = "barrnap")
    tempfile_rnafa = paste0(tempfilebase, ".", kingdom, ".fa")
    tempfile_rnagff = paste0(tempfilebase, ".", kingdom, ".gff")
    tempfile_16S = paste0(tempfilebase, ".", kingdom, ".16S.txt")

    system2("barrnap",
            args = c("--quiet",
                     "--threads", threads,
                     "--kingdom", kingdom,
                     "--outseq", tempfile_rnafa
            ),
            stdin = genome_file,
            stdout = tempfile_rnagff,
            stderr = "/dev/null"
    )
    system2("grep",
            args = c("-e", "Name=16S_rRNA"),
            stdin = tempfile_rnagff,
            stdout = tempfile_16S,
            stderr = "/dev/null"
    )
    #file.remove(tmp.file)
    result = list(rRNA_fa = tempfile_rnafa,
                  rRNA_gff = tempfile_rnagff,
                  tempfile_16S = tempfile_16S)
    return(result)
  }
}

run.bedtools <- function(genome_file = system.file("extdata", "2503283023.fna", package = "microtrait", mustWork = TRUE),
                         bed16S_file) {
  # Check arguments
  genome_file = tryCatch({
    assertthat::assert_that(file.exists(genome_file))
    genome_file
  },
  error = function(e) {
    message("Genome file ", `genome_file`, " doesn't exist.")
    print(e)
  }
  )

  bed16S_file = tryCatch({
    assertthat::assert_that(file.exists(bed16S_file))
    bed16S_file
  },
  error = function(e) {
    message("16S bed file ", `bed16S_file`, " doesn't exist.")
    print(e)
  }
  )
  #checkmate::qexpect(genome_file, "S1", info = "info", label = "Genome fasta file")
  #checkmate::qassert(faa_file, "S1")
  #checkmate::assertChoice(mode, c("single", "meta"))
  #checkmate::qassert(maskN, "B")
  #checkmate::qassert(bypassSD, "B")
  #n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", parallel::detectCores(), "]")))
  #prodigal -i ./2775507255/2775507255.fna -o ./2775507255.prodigal -a ./2775507255.prodigal.faa -p single -q -f gff -m -n -g 11
  result <- NULL

  #bedtools getfasta -s -name -fi genome_file -bed tempfile_16S
  # -fo ./output/genomes/desulfotomaculum_kuznetsovii/Desulfotomaculum_kuznetsovii_dsm_6115_ASM21470v1_dna_toplevel/barrnap/barrnap_results_assigned.fa

  if(available.external("bedtools")){
    message("Running bedtools for ", fs::path_file(genome_file), " and ", fs::path_file(bed16S_file))
    # to do: reimplement with processx
    tempfilebase = tempfile(pattern = "barrnap")
    tempfile_16Sfa = paste0(tempfilebase, ".16S.fa")

    system2("bedtools",
            args = c("getfasta",
                     "-s -name",
                     "-fi", genome_file,
                     "-bed", bed16S_file,
                     "-fo", tempfile_16Sfa
            ),
            stdout = "/dev/null",
            stderr = "/dev/null"
    )
    #file.remove(tmp.file)
    return(tempfile_16Sfa)
  }
}

available.external <- function(what) {
  if(what == "hmmsearch"){
    chr <- NULL
    try(chr <- system2("hmmsearch", args = "-h", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('hmmsearch was not found by R.',
                 'Please install hmmsearch from: http://hmmer.org/download.html',
                 'After installation, make sure hmmsearch can be run from a shell/terminal, ',
                 'using the command \'hmmsearch --help\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  if(what == "hmmfetch"){
    chr <- NULL
    try(chr <- system2("hmmfetch", args = "-h", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('hmmfetch was not found by R.',
                 'Please install hmmfetch from: http://hmmer.org/download.html',
                 'After installation, make sure hmmfetch can be run from a shell/terminal, ',
                 'using the command \'hmmfetch --help\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  if(what == "hmmpress") {
    chr <- NULL
    try(chr <- system2("hmmpress", args = "-h", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('hmmpress was not found by R.',
                 'Please install hmmpress from: http://hmmer.org/download.html',
                 'After installation, make sure hmmpress can be run from a shell/terminal, ',
                 'using the command \'hmmpress --help\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  if(what == "prodigal") {
    chr <- NULL
    try(chr <- system2("prodigal", args = "-v", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('prodigal was not found by R.',
                 'Please install barrnap from: https://github.com/hyattpd/Prodigal',
                 'After installation, make sure prodigal can be run from a shell/terminal, ',
                 'using the command \'prodigal -h\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  if(what == "tRNAscan") {
    chr <- NULL
    try(chr <- system2("tRNAscan-SE", args = "-h", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('tRNAscan-SE was not found by R.',
                 'Please install tRNAscan-SE from: http://lowelab.ucsc.edu/tRNAscan-SE/',
                 'After installation, make sure tRNAscan-SE can be run from a shell/terminal, ',
                 'using the command \'tRNAscan-SE -h\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  if(what == "cmsearch") {
    chr <- NULL
    try(chr <- system2("cmsearch", args = "-h", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('cmsearch was not found by R.',
                 'Please install cmsearch from: http://eddylab.org/infernal/',
                 'After installation, make sure cmsearch can be run from a shell/terminal, ',
                 'using the command \'cmsearch -h\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  if(what == "barrnap") {
    chr <- NULL
    try(chr <- system2("barrnap", args = "-v", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('barrnap was not found by R.',
                 'Please install barrnap from: https://github.com/tseemann/barrnap',
                 'After installation, make sure barrnap can be run from a shell/terminal, ',
                 'using the command \'barrnap -h\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  if(what == "bedtools") {
    chr <- NULL
    try(chr <- system2("bedtools", args = "--version", stdout = TRUE, stderr= TRUE), silent = TRUE)
    if(is.null(chr)){
      stop(paste('bedtools was not found by R.',
                 'Please install bedtools from: https://github.com/arq5x/bedtools2/releases',
                 'After installation, make sure bedtools can be run from a shell/terminal, ',
                 'using the command \'bedtools --help\', then restart the R session.', sep = '\n'))
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}

countseq.fasta <- function(fastafile) {
  nseq = as.numeric(system(paste0("grep \">\" ", fastafile, "|wc -l"), intern = T))
  nseq
}
