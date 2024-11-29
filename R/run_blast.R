#' @title Run BLAST
#'
#' @description Run BLAST for a single sequence and return raw results.
#'   For formatted results, use the function *`BLASTr::get_blastn_results()`*
#'
#' @param asv Vector of sequences to be BLASTed.
#' @param db_path Complete path do formatted BLAST database.
#' @param num_threads Number of threads to run BLAST on.
#'   Passed on to BLAST+ argument `-num_threads`.
#' @param perc_id Lowest identity percentage cutoff.
#'   Passed on to BLAST+ `-perc_identity`.
#' @param perc_qcov_hsp Lowest query coverage per HSP percentage cutoff.
#' Passed on to BLAST+ `-qcov_hsp_perc`.
#' @param num_alignments Number of alignments to retrieve from BLAST. Max = 6.
#' @param blast_type One of the available BLAST+ search engines,
#' #'   one of: `c("blastn", "blastp", "blastx", "tblastn", "tblastx")`.
#' @param verbose Should condathis::run() internal command be shown?
#' @param env_name The name of the conda environment with the parameter (i.e. "blast-env")
#' @return Unformatted BLAST results.
#'   For results formatted as tibble, please use `BLASTr::get_blast_results()`
#'
#' @examples
#' \dontrun{
#' blast_res <- run_blast(
#'   asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
#'     AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
#'   db_path = fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta"),
#'   perc_id = 80,
#'   num_thread = 1,
#'   perc_qcov_hsp = 80,
#'   num_alignments = 2
#' )
#' }
#'
#' @export
run_blast <- function(asv,
                      db_path,
                      num_alignments = 4L,
                      num_threads = 1L,
                      blast_type = "blastn",
                      # gapopen = 5,
                      # gapextend = 2,
                      # task = "blastn",
                      perc_id = 80L,
                      perc_qcov_hsp = 80L,
                      verbose = FALSE,
                      env_name = "blast-env") {
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
  # if (is.null(num_threads)) {
  #   num_threads <- getOption(
  #     "BLASTr.num_threads",
  #     default = 1
  #   )
  # }

  rlang::check_required(asv)
  rlang::check_required(db_path)
  check_cmd(blast_type, env_name = env_name, verbose = verbose)

  query_path <- fs::file_temp("blast_input_", ext = "fasta")
  base::cat(asv, file = query_path)

  blast_res <- condathis::run_bin(
    blast_type,
    "-db", db_path,
    "-query", query_path,
    "-outfmt", "6 std qcovhsp staxid",
    # "-gapopen", gapopen,
    # "-gapextend", gapextend,
    # "-task", task,
    "-max_hsps", "1",
    "-perc_identity", perc_id,
    "-qcov_hsp_perc", perc_qcov_hsp,
    "-num_threads", as.character(num_threads),
    "-num_alignments", as.character(num_alignments),
    env_name = env_name,
    verbose = verbose,
    error = "continue"
  )

  if (isTRUE(blast_res$status != 0L)) {
    error_msg_list <- stringr::str_extract_all(blast_res$stderr, "Error:.*")

    if (length(unlist(error_msg_list)) > 0) {
      error_msg_vector <- unlist(error_msg_list)
      names(error_msg_vector) <- c(rep("x", times = length(error_msg_vector)))
    } else {
      error_msg_vector <- blast_res$stderr
    }

    cli::cli_abort(
      message = error_msg_vector,
      class = "blastr_error_blast_run"
    )
  }

  return(blast_res)
}
