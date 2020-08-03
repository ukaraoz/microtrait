#' @examples
#' rds.file = "/Users/ukaraoz/Work/microtrait/code/microtrait-analysis/ncbimags/rds/PRJNA288027/GCA_001792575.1_ASM179257v1_genomic.fna.microtrait.rds"
#' fetch_microtrait_results(rds.file, id = "GCA_001792575.1_ASM179257v1", type = "trait")
fetch_microtrait_results <- function(rds.file, id, type = "trait") {
  temp = readRDS(rds.file)
  if(type == "trait"){
    result = temp$all_traits %>% dplyr::select(c("trait-display-short", "n"))
  }
  if(type == "gene"){
    result = temp$genes_detected %>% tbl_df
  }
  if(type == "rule"){
    result = temp$rules_asserted
  }
  if(type == "time"){
    result = data.frame(time = as.numeric(gsub(".*: | sec.*", "", temp$time_log[[5]])))
  }
  result = result %>%
    tibble::add_column(id = id) %>%
    dplyr::select(id, everything()) %>% as.data.frame
  return(result)
}

