% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extern.R
\name{run.hmmsearch}
\alias{run.hmmsearch}
\title{Run hmmsearch (HMMER 3)}
\usage{
run.hmmsearch(
  faa_file = system.file("extdata", "2619619645.genes.faa", package = "microtrait",
    mustWork = TRUE),
  hmm = NULL,
  domtblout_file = NULL,
  n_threads = 1
)
}
\arguments{
\item{hmm}{A hmm file. Must be pressed (see hmmpress from HMMER manual).}

\item{domtblout_file}{}

\item{n_threads}{An \code{integer}. The number of cores to use.}

\item{fasta}{A protein fasta file.}
}
\value{
The path to a temporary file where the hmmsearch output is placed.

hmmsearch --cpu $ncpus --domtblout $faafile".microtrait.domtblout" $hmm $faafile > /dev/null
}
\description{
Takes a fasta file and a Hidden Markov Model profile and
performs a search of the former over the latter.
}
