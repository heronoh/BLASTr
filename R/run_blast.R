#' @title Run BLAST
#'
#' @description Run BLAST for a single sequence and return raw results. For formated results, use the function *BLASTr::getblastn_results()*
#'
#' @param asvs Vector of sequences to be BLASTed.
#' @param db_path Complete path do formatted BLAST database.
#' @param num_threads Number of threads to run BLAST on. Passed on to BLAST+ argument _-num_threads_.
#' @param perc_ID Lowest identity percentage cutoff. Passed on to BLAST+ _-perc_identity_.
#' @param perc_qcov_hsp Lowest query coverage per HSP percentage cutoff. Passed on to BLAST+ _-qcov_hsp_perc_.
#' @param num_alignments Number of alignments to retrieve from BLAST. Max = 6.
#' @param blast_type One of the available BLAST+ search engines: "blastn", "blastp", "blastx", "tblastn", "tblastx".
#'
#' @return Unformatted BLAST results. For results formatted as tibble, please use *BLASTr::get_blast_results()*
#'
#' @examples
#' blast_res <- BLASTr::run_blast(asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
#'                                 db_path = "/data/databases/nt/nt",
#'                                 perc_ID = 80,
#'                                 num_thread = 1,
#'                                 perc_qcov_hsp = 80,
#'                                 num_alignments = 2)
#'
#' @export

run_blast <- function(
  asv,
  db_path,
  num_alignments,
  num_threads,
  blast_type,
  perc_ID,
  perc_qcov_hsp
  ) {
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

  # blast_bin <- check_bin("blastn")

  # blast_type <- "blastn"

  blast_bin <- check_bin(blast_type)

  blast_cmd <- "{blast_bin} -db {db_path} -outfmt '6 std qcovhsp' -max_hsps 1 -perc_identity {perc_ID} -qcov_hsp_perc {perc_qcov_hsp} -num_threads {as.character(num_threads)} -num_alignments {as.character(num_alignments)}"

  blast_cmd_in <- paste0("echo -e '>seq1\n{asv}' | ", blast_cmd)

  # print(blast_cmd_in)

  # previous version of completeBLAST command inside the function
  # blast_cmd_original <- "echo -e '>seq1\n{asv}' | blastn -db {db_path} -outfmt '6 std qcovhsp' -max_hsps 1 -perc_identity 80 -qcov_hsp_perc 80 -num_threads {as.character(num_thread)} -num_alignments {as.character(num_alignments)}"

  blast_res <- shell_exec(
    cmd = blast_cmd_in
  )
  return(blast_res)
}
