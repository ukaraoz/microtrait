readFASTA <- function (file = system.file("protseq/P00750.fasta", package = "protr"),
                       legacy.mode = TRUE, seqonly = FALSE) {
  # Read the FASTA file as a vector of strings
  lines <- readLines(file)

  # Remove comment lines starting with a semicolon ';'
  if (legacy.mode) {
    comments <- grep("^;", lines)
    if (length(comments) > 0) {
      lines <- lines[-comments]
    }
  }

  # Get the line numbers where sequences names are
  ind <- which(substr(lines, 1L, 1L) == ">")

  # Compute the total number of sequences
  nseq <- length(ind)

  if (nseq == 0) stop("no line starting with a > character found")

  # Localize sequence data
  start <- ind + 1
  end <- ind - 1
  end <- c(end[-1], length(lines))

  # Read sequences
  sequences <- lapply(
    seq_len(nseq),
    function(i) paste(lines[start[i]:end[i]], collapse = "")
  )

  if (seqonly) return(sequences)

  # Read sequence names
  nomseq <- lapply(seq_len(nseq), function(i) {
    firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
    substr(firstword, 2, nchar(firstword))
  })

  # Give the sequences names to the list elements
  names(sequences) <- nomseq
  sequences
}
