#' Prepare dbcan database (download and subselect)
#'
#' @import futile.logger
download.dbcan <- function(dbcan_version = 8, dbcanhmmdb_selectids_file, dbcanhmmdb_file) {
  futile.logger::flog.info("Downloading dbcan hmm database")

  dbcan_hmmdb_url = paste("http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V", dbcan_version, ".txt", sep = "")
  downloaded_file <- file.path(tempdir(), paste0("dbcan.v", dbcan_version, ".txt"))

  download.file(dbcan_hmmdb_url, downloaded_file)

  futile.logger::flog.info("Subsetting dbcan hmm database")
  if(available.external("hmmfetch")) {
    system(paste("hmmfetch",
                 "-f", downloaded_file,
                 dbcanhmmdb_selectids_file,">", dbcanhmmdb_file
    )
    )
  }
}

#' Prepare microtrait database (download and subselect)
#' @import futile.logger piggyback fs
#' @importFrom R.utils gunzip
download.microtrait <- function(microtraithmmdb_file) {
  futile.logger::flog.info("Downloading microtrait hmm database")

  #piggyback::pb_download("microtrait.hmmdb.gz",
  #                       repo = "ukaraoz/microtrait-hmmtest",
  #                       #repo = "ukaraoz/test",
  #                       dest = fs::path_dir(microtraithmmdb_file))

  microtrait_hmmdb_url = "https://github.com/ukaraoz/microtrait-hmmtest/releases/download/latest/microtrait.hmmdb.gz"
  download.file(microtrait_hmmdb_url,
                destfile = microtraithmmdb_file)

  microtraithmmdb_unzippedfile = R.utils::gunzip(microtraithmmdb_file, remove = F, overwrite = T)
  return(microtraithmmdb_unzippedfile[1])
}

#' Prepare microtrait database (download and subselect)
#' @import futile.logger
#' @export prep.hmmmodels
prep.hmmmodels <- function(output_dir=system.file("extdata", package = "microtrait")) {
  hmmpress.extensions = c("h3f", "h3i", "h3m", "h3p")
  hmmpress_dir = file.path(output_dir, "hmm", "hmmpress")

  # start clean
  unlink(file.path(output_dir, "hmm", "dbcan", "*"))
  unlink(file.path(output_dir, "hmm", "microtrait-hmmdb", "*"))
  unlink(file.path(output_dir, "hmm", "hmmpress", "*"))

  dbcan_version = 8
  dbcanhmmdb_selectids_file = file.path(output_dir, "dbcan.selectids.txt")
  dbcanhmmdb_file = file.path(output_dir, "hmm", "dbcan", paste0("dbcan.select.v", dbcan_version, ".hmmdb"))
  microtraithmmdb_file = file.path(output_dir, "hmm", "microtrait-hmmdb", "microtrait.hmmdb.gz")

  #if(!file.exists(dbcanhmmdb_file)) {
  download.dbcan(dbcanhmmdb_selectids_file = dbcanhmmdb_selectids_file,
                 dbcanhmmdb_file = dbcanhmmdb_file)
  futile.logger::flog.info("Preparing dbcan hmm database")
  if(available.external("hmmpress")) {
    system(paste("hmmpress", dbcanhmmdb_file))
  }
  file.copy(paste(dbcanhmmdb_file, hmmpress.extensions, sep = "."), hmmpress_dir)
  file.remove(paste(dbcanhmmdb_file, hmmpress.extensions, sep = "."))

  # microtraithmmdb_file = "/Users/ukaraoz/Work/microtrait/code/microtrait/inst/extdata/hmm/microtrait-hmmdb/microtrait.hmmdb.tar.gz"
  microtraithmmdb_unzippedfile = download.microtrait(microtraithmmdb_file)
  futile.logger::flog.info("Preparing microtrait hmmdb")
  system(paste("hmmpress", microtraithmmdb_unzippedfile))

  file.copy(paste(microtraithmmdb_unzippedfile, hmmpress.extensions, sep = "."), hmmpress_dir)
  file.remove(paste(microtraithmmdb_unzippedfile, hmmpress.extensions, sep = "."))
}
