run.pullseq <- function(faa_file,
                        names,
                        out_file,
                        regex = F, 
                        pullseq = "/opt/pullseq/src/pullseq") {
  message(paste0("Writing to ", out_file))
  if(regex == F) { system2(pullseq,
          args = c("--input", faa_file,
                   "--names", names
          ),
          stdout = out_file,
          stderr = "/dev/null"
         )
  }else if(regex == T) {
    system2(pullseq,
          args = c("--input", faa_file,
                   "--regex", names
          ),
          stdout = out_file,
          stderr = "/dev/null"
         )
  }
}

detect.loci <- function(faa_file,
                        hmm_file,
                        tblout_file,
                        usearch_db,
                        usearch_outfile,
                        option = c("hmmscan", "usearch", "both")) {
  if(!missing(option) & length(option)>1) stop("Only one 'method' allowed.")
  option <- match.arg(option)
  if(option=="hmmscan") {
    tblout = run.hmmsearch(faa_file = faa_file, 
                               tblout_file = tblout_file,
                               hmm_file = hmm_file)
    result = tblout %>% select(gene_name) %>% distinct %>% tibble::add_column(method = "hmmscan")
  } else if(option=="usearch") {
    usearch_out = run.usearch(faa_file = faa_file, 
                              usearch_db = usearch_db,
                              usearch_outfile = usearch_outfile)
    result = usearch_out %>% select(gene_name) %>% distinct %>% tibble::add_column(method = "usearch")
    #return(cat(faa_file, "\t", "usearch", "\n"))
  } else {
    message("Running hmmsearch")
    tblout = run.hmmsearch(faa_file = faa_file, 
                               tblout_file = tblout_file,
                               hmm_file = hmm_file)
    message("Running usearch")
    usearch_out = run.usearch(faa_file = faa_file, 
                              usearch_db = usearch_db,
                              usearch_outfile = usearch_outfile)
    result1 = tblout %>% select(gene_name) %>% distinct %>% tibble::add_column(method = "hmmscan")
    result2 = usearch_out %>% select(gene_name) %>% distinct %>% tibble::add_column(method = "usearch")
    result = result1 %>% tibble::add_row(result2)
  }
  return(result)
}

scan.faa <- function(faa_file,
                     hmm_files = hmm.files,
                     out_dir = dirname(file.path(faa_file))) {
  print(paste0("Scanning ", faa_file))
  results = data.frame(gene_name = character(),
                     method = character(),
                     faa = character(),
                     loci = character()
                     ) %>% tbl_df
  for(i in 1:length(hmm_files)) {
    loci = sub(".hmm", "", basename(hmm_files[i]))
    detected = detect.loci(faa_file = faa_file, 
                           hmm_file = hmm_files[i],
                           tblout_file = file.path(out_dir, paste0(basename(faa_file), ".", loci, ".tblout")),
                           option = "hmmscan")
    detected = detected %>% 
                tibble::add_column(faa = faa_file) %>% 
                tibble::add_column(loci = sub("\\..*", "", basename(hmm_files[i])))
    results = results %>%
                tibble::add_row(detected)
  }
  results
}

run.usearch <- function(faa_file = NULL,
                        usearch_db = usearch_db,
                        usearch_outfile = NULL,
                        evalue = 0.01,
                        cov_target_threshold = 70) {
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
  #message("Running usearch for ", fs::path_file(faa_file))
  #message(command)
  #command = paste("hmmsearch", "--domtblout", domtblout_file, "--cpu", n_threads,"--noali",
  #                 hmm_file, faa_file, sep = " ")

  system2("/opt/usearch",
          args = c("-ublast", faa_file,
                   "-db", usearch_db,
                   "-evalue", evalue,
                   "-threads", "6",
                   "-userout", usearch_outfile,
                   "-userfields", "query+target+id+alnlen+ql+tl+mism+opens+qlo+qhi+tlo+thi+evalue+bits"),
          stdout = "/dev/null",
          stderr = "/dev/null"
  )

  col_types <- readr::cols(gene_name = readr::col_character(),
                           target = readr::col_character(),
                           id = readr::col_double(),
                           alnlen = readr::col_integer(),
                           ql = readr::col_integer(),
                           tl = readr::col_integer(),
                           mism = readr::col_integer(),
                           opens = readr::col_integer(),
                           qlo = readr::col_integer(),
                           qhi = readr::col_integer(),
                           tlo = readr::col_integer(),
                           thi = readr::col_integer(),
                           evalue = readr::col_double(),
                           bits = readr::col_integer()
  )
  usearch_out = readr::read_lines(usearch_outfile)
  if(length(usearch_out) != 0) {
    usearch_out = usearch_out %>%
    readr::read_tsv(col_names=names(col_types$cols), na='-') %>%
    readr::type_convert(col_types=col_types) %>%
    dplyr::group_by(gene_name,target) %>%
    dplyr::arrange( desc(bits) ) %>%
    # Pick the top 1 value
    dplyr::slice(1) %>%
    dplyr::mutate(cov_query = ((qhi-qlo)/ql)*100,
                  cov_target = ((thi-tlo)/tl)*100) %>%
    dplyr::select(gene_name, target, id, cov_query, cov_target, evalue, bits) %>%
    dplyr::filter(cov_target >= cov_target_threshold) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene_name = sub("(.?) #.*", "\\1", gene_name))
  }
  #file.remove(tmp.file)
  return(usearch_out)
}


run.hmmsearch <- function(faa_file = NULL,
                          hmm_file = NULL,
                          tblout_file = NULL,
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
  
  #hmmsearch --cpu $ncpus --domtblout $faafile".microtrait.domtblout" $hmm $faafile > /dev/null
  #message("Running hmmsearch for ", fs::path_file(faa_file))
  #message(command)
  #command = paste("hmmsearch", "--domtblout", domtblout_file, "--cpu", n_threads,"--noali",
  #                 hmm_file, faa_file, sep = " ")

  system2("hmmsearch",
          args = c("--tblout", tblout_file,
                   "--cpu", n_threads,
                   "--cut_tc --noali",
                   hmm_file, faa_file
          ),
          stdout = "/dev/null",
          stderr = "/dev/null"
  )
  #file.remove(tmp.file)
  tblout = read.tblout(tblout_file)

  return(tblout)
}

read.phylodist <- function(phylodist_file) {
  print(paste0("Reading ", phylodist_file))
  col_types <- readr::cols(
      gene_id = readr::col_character(),
      homolog_gene_oid = readr::col_character(),
      homolog_taxon_oid = readr::col_character(),
      percent_identity= readr::col_double(),
      lineage = readr::col_character())
  phylodist = read_tsv(phylodist_file, 
                  col_names = c("gene_id", "homolog_gene_oid", "homolog_taxon_oid", "percent_identity", "lineage"), 
                  col_types = col_types)
  return(phylodist)
}

read.ggkbasephylodist <- function(phylodist_file) {
  print(paste0("Reading ", phylodist_file))
  col_types <- readr::cols(
      `Contig name` = readr::col_character(),
      `Size (bp)` = readr::col_double(),
      Coverage = readr::col_double(),
      `GC %` = readr::col_double(),
      `Taxonomy winner` = readr::col_character(),
      `Winner %` = readr::col_double(),
      `Species winner` = readr::col_character(),
      `Species winner %` = readr::col_double(),
      `Genus winner` = readr::col_character(),
      `Genus winner %` = readr::col_double(),
      `Order winner` = readr::col_character(),
      `Order winner %` = readr::col_double(),
      `Class winner` = readr::col_character(),
      `Class winner %` = readr::col_double(),
      `Phylum winner` = readr::col_character(),
      `Phylum winner %` = readr::col_double(),
      `Domain winner` = readr::col_character(),
      `Domain winner %` = readr::col_double())
  col_names = c("contig_name", "size", "coverage", "gc", "taxonomy_winner", "winner_percent", 
                                "species_winner", "species_winner_percent", "genus_winner", "genus_winner_percent", "order_winner", "order_winner_percent",
                                "class_winner", "class_winner_percent", "phylum_winner", "phylum_winner_percent", "domain_winner", "domain_winner_percent")
  phylodist = read_tsv(phylodist_file, col_types = col_types)
  colnames(phylodist) = col_names
  return(phylodist)
}
