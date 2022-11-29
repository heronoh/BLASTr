#1b - run_blast ----
###                  function to run blast for each ASV/OTU                  ###
#' Run BLAST
#'
#' Run BLAST for a single sequence
#'
#' @param asv Sequence to be BLASTed
#' @param db_path Path to BLAST database
#' @param num_alignments Number of alignments to be returned.
#'   Default: the 3 best alignments.
#' @param num_thread Number of threads to be used
#' @param perc_ID Percent identity to be used
#' @param perc_qcov_hsp Percent query coverage to be used
#'
#' @export
run_blast <- function(
  asv,
  db_path,
  num_alignments = 3,
  num_thread,
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
  if (is.null(num_thread)) {
    num_thread <- getOption(
      "BLASTr.num_thread",
      default = 1
    )
  }

  blastn_bin <- check_bin("blastn")

  blast_cmd <- "{blastn_bin} -db {db_path} -outfmt '6 std qcovhsp' -max_hsps 1 -perc_identity {perc_ID} -qcov_hsp_perc {perc_qcov_hsp} -num_threads {as.character(num_thread)} -num_alignments {as.character(num_alignments)}"

  blast_cmd_in <- paste0("echo -e '>seq1\n{asv}' | ", blast_cmd)

  # print(blast_cmd_in)

  #previous version of completeBLAST command inside the function
  # blast_cmd_original <- "echo -e '>seq1\n{asv}' | blastn -db {db_path} -outfmt '6 std qcovhsp' -max_hsps 1 -perc_identity 80 -qcov_hsp_perc 80 -num_threads {as.character(num_thread)} -num_alignments {as.character(num_alignments)}"

  blast_res <- shell_exec(
    cmd = blast_cmd_in
  )
  return(blast_res)
}
