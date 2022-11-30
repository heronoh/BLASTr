#2 - Running BLASTr function with supplied ASVs ----
#2a - Set parallelization ----


#' Run BLAST on parallel
#'
#' Run parallel BLAST for set of sequences
#'
#' @param out_file Complete path to output file.
#' @param out_RDS Complete path to output RDS file.
#' @param total_cores Total available cores to parallelize BLAST. Chech your max with _future::availableCores()_
#' @param blast_type BLAST+ executable to be used on search
#'
#' @inheritParams run_blast
#'
#' @return A tibble with the BLAST tabular output.
#'
#' @export
#'
parallel_blast <- function(
  asvs,
  db_path,
  out_file,
  out_RDS,
  num_thread,
  # blast_type,
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
                                       num_thread = 1,
                                       # blast_type = blast_type,
                                       num_alignments = num_alignments,
                                       db_path = db_path,
                                       perc_ID = perc_ID,
                                       perc_qcov_hsp = perc_qcov_hsp
                                       )
    } else {
      blast_res <- purrr::map_dfr(.x = asvs,
                                  .f = get_blast_results,
                                  num_thread = 1,
                                  # blast_type = blast_type,
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

