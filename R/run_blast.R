#' @title Run BLAST
#'
#' @description Run BLAST for a single sequence and return raw results.
#'   For formated results, use the function *BLASTr::getblastn_results()*
#'
#' @param asvs Vector of sequences to be BLASTed.
#' @param db_path Complete path do formatted BLAST database.
#' @param num_threads Number of threads to run BLAST on.
#'   Passed on to BLAST+ argument _-num_threads_.
#' @param perc_id Lowest identity percentage cutoff.
#'   Passed on to BLAST+ _-perc_identity_.
#' @param perc_qcov_hsp Lowest query coverage per HSP percentage cutoff.
#' Passed on to BLAST+ _-qcov_hsp_perc_.
#' @param num_alignments Number of alignments to retrieve from BLAST. Max = 6.
#' @param blast_type One of the available BLAST+ search engines,
#'   one of: "blastn", "blastp", "blastx", "tblastn", "tblastx".
#'
#' @return Unformatted BLAST results.
#'   For results formatted as tibble, please use `BLASTr::get_blast_results()`
#'
#' @examples
#' blast_res <- BLASTr::run_blast(
#'   asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
#'     AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
#'   db_path = "/data/databases/nt/nt",
#'   perc_id = 80,
#'   num_thread = 1,
#'   perc_qcov_hsp = 80,
#'   num_alignments = 2
#' )
#'
#' @export
run_blast <- function(asv,
                      db_path,
                      num_alignments = 4,
                      num_threads = 1,
                      blast_type = "blastn",
                      perc_id = 80,
                      perc_qcov_hsp = 80) {
  #   if (is.null(db_path)) {
  #   db_path <- getOption(
  #     "BLASTr.db_path",
  #     default = NULL
  #   )
  # }
  # if (is.null(db_path)) {
  #   cli::cli_abort(
  #     message = "No BLAST database provided."
  #   )
  # }
  if (is.null(num_threads)) {
    num_threads <- getOption(
      "BLASTr.num_threads",
      default = 1
    )
  }
  blast_bin <- BLASTr:::check_bin(blast_type)
  # rlang::inform(blast_bin)
  blast_cmd <- "{blast_bin} -db {db_path} -outfmt '6 std qcovhsp staxid' -max_hsps 1 -perc_identity {perc_id} -qcov_hsp_perc {perc_qcov_hsp} -num_threads {as.character(num_threads)} -num_alignments {as.character(num_alignments)}"
  blast_cmd_in <- paste0("echo -e '>seq1\n{asv}' | ", blast_cmd)
  blast_res <- BLASTr:::shell_exec(
    cmd = blast_cmd_in
  )
  return(blast_res)
}
