#' microtraithmm_fromrules.
#'
#' HMMs (Hidden Markov Models) underlying microTrait pipeline.
#' This data provides information on the database cross references for
#' each HMM and the performance of the HMM using KEGG orthologs as a
#' benchmark.
#'
#' @format A data frame with 13 variables:
#' \describe{
#'   \item{microtraithmm-name}{}
#'   \item{microtraithmm-dbxref_kegg}{}
#'   \item{microtraithmm-dbxref_ec}{}
#'   \item{microtraithmm-dbxref_TC}{}
#'   \item{microtraithmm-description}{}
#'   \item{microtraithmm-notes}{}
#'   \item{score}{}
#'   \item{fscore}{}
#'   \item{tpr}{}
#'   \item{tnr}{}
#'   \item{fpr}{}
#'   \item{fnr}{}
#'   \item{accuracy}{}
#' }
"hmms_fromrules"

#' microtraitrule_table.
#'
#' Rules underlying microTrait pipeline.
#' Each rule is a logical rule using presence/absence of microTrait HMMs.
#'
#' @format A data frame with 13 variables:
#' \describe{
#'   \item{microtraithmm-name}{}
#'   \item{microtraithmm-rule}{}
#'   \item{microtraithmm-name-quoted}{}
#'   \item{microtraithmm-ruleunwrapped}{}
#' }
"rules"
