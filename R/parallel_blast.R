#' @title Run BLAST on parallel
#'
#' @description Run parallel BLAST for set of sequences
#'
#' @inheritParams run_blast
#'
#' @param out_file Complete path to output file
#' @param out_RDS Complete path to output RDS file.
#' @param total_cores Total available cores to parallelize BLAST. Chech your max with *future::availableCores()*
#' @param blast_type BLAST+ executable to be used on search.
#'
#' @return A tibble with the BLAST tabular output.
#'
#' @examples
#' blast_res <- BLASTr::parallel_blast(asvs = ASVs_test[1],
#'                                  db_path = "/data/databases/nt/nt",
#'                                  out_file = "~/blast_out.csv",
#'                                  out_RDS = "~/blast_out.RDS",
#'                                  total_cores = 77,
#'                                  perc_ID = 80,
#'                                  num_thread = 1,
#'                                  perc_qcov_hsp = 80,
#'                                  num_alignments = 4)
#'
#' @export

parallel_blast <- function(
  asvs,
  db_path,
  out_file,
  out_RDS,
  num_threads,
  blast_type,
  total_cores,
  perc_ID,
  perc_qcov_hsp,
  num_alignments
  ){

  # TODO: Convert ASVs to vector, if needed


  # if (is.null(db_path)) {
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
  # if (is.null(num_thread)) {
  #   num_thread <- getOption(
  #     "BLASTr.num_thread",
  #     default = 1
  #   )
  # }
  if (total_cores > 1) {
    future::plan(future::multisession(),
                 workers = total_cores)

    blast_res <- furrr::future_map_dfr(.x = asvs,
                                       .f = get_blast_results,
                                       .options = furrr::furrr_options(seed = NULL),
                                       num_threads = 1,
                                       blast_type = blast_type,
                                       num_alignments = num_alignments,
                                       db_path = db_path,
                                       perc_ID = perc_ID,
                                       perc_qcov_hsp = perc_qcov_hsp
                                       )
    } else {
      blast_res <- purrr::map_dfr(.x = asvs,
                                  .f = get_blast_results,
                                  num_threads = 1,
                                  blast_type = blast_type,
                                  num_alignments = num_alignments,
                                  db_path = db_path,
                                  perc_ID = perc_ID,
                                  perc_qcov_hsp = perc_qcov_hsp
                                )

  }
  if (!is.na(out_file)) {
    readr::write_csv(
      x = blast_res,
      file = out_file,
      append = FALSE
    )
  }
  if (!is.na(out_RDS)) {
    readr::write_rds(
      x = blast_res,
      file = out_RDS
    )
  }
  return(blast_res)
}
