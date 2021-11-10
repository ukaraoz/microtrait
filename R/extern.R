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

  hmm = match.arg(hmm, c("microtrait", "dbcan"))
  hmm_file = switch(hmm,
                    #"microtrait" = file.path("/Users/ukaraoz/Work/microtrait/code/microtrait", "inst/extdata/hmm/hmmpress/microtrait.hmmdb"),
                    #"dbcan" = file.path("/Users/ukaraoz/Work/microtrait/code/microtrait", "inst/extdata/hmm/hmmpress/dbcan.select.v8.hmmdb")
                    "microtrait" = file.path(system.file(package="microtrait"), "extdata/hmm/hmmpress/microtrait.hmmdb"),
                    "dbcan" = file.path(system.file(package="microtrait"), "extdata/hmm/hmmpress/dbcan.select.v8.hmmdb")
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
    message("Running hmmsearch for ", fs::path_file(faa_file))
    #message(command)
    if(hmm == "microtrait") {
      system2("hmmsearch",
              args = c("--domtblout", domtblout_file,
                       "--cpu", n_threads,
                       "--noali",
                       "--cut_tc ", # CRUCIAL, otherwise doesn't use model specific thresholds
                       hmm_file, faa_file
              ), stdout = "/dev/null")
    } else if(hmm == "dbcan") {
      system2("hmmsearch",
              args = c("--domtblout", domtblout_file,
                       "--cpu", n_threads,
                       "--noali",
                       hmm_file, faa_file
              ), stdout = "/dev/null")
    }

    #file.remove(tmp.file)
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
#' @return
#'
#' @note The prodigal software must be installed on the system for this function to work, i.e. the command
#' \samp{system("prodigal -h")} must be recognized as a valid command if you run it in the Console window.
#'
#' @seealso
run.prodigal <- function(genome_file = system.file("extdata/examples/2619619645/in", "2619619645.genes.fna", package = "microtrait", mustWork = TRUE),
                         faa_file,
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

  result <- NULL
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
                     "-a", faa_file,
                     "-o", tempfile(pattern = "prodigal", fileext = "gff")
            ),
            stdout = "/dev/null",stderr = "/dev/null"
    )
    #file.remove(tmp.file)
    return(faa_file)
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

  if(what == "hmmpress"){
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

  if(what == "prodigal"){
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
}

countseq.fasta <- function(fastafile) {
  nseq = as.numeric(system(paste0("grep \">\" ", fastafile, "|wc -l"), intern = T))
  nseq
}
