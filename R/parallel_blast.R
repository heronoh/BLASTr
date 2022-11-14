#2 - Running BLASTr function with supplied ASVs ----
#2a - Set parallelization ----


#' Run Parallel BLAST
#'
#' Run parallel BLAST for set of sequences
#'
#' @param asvs Vector of sequences to be BLASTed
#' @param out_file Path to output file
#' @param out_rds Path to output RDS file
#'
#' @inheritParams run_blast
#'
#' @export
parallel_blast <- function(
  asvs,
  db_path,
  out_file,
  out_rds,
  num_thread,
  perc_id,
  perc_qcov_hsp
){
  # TODO: Convert ASVs to vector, if needed
  if (is.null(db_path)) {
    db_path <- getOption(
      "BLASTr.db_path",
      default = NULL
    )
  }
  if (is.null(db_path)) {
    cli::cli_abort(
      message = "No BLAST database provided."
    )
  }
  if (is.null(num_thread)) {
    num_thread <- getOption(
      "BLASTr.num_thread",
      default = 1
    )
  }
  if (num_thread > 1) {
    future::plan(future::multisession(), workers = num_thread)
    blast_res <- furrr::future_map_dfr(
      asvs,
      get_blast_results,
      num_thread = 1,
      .options = furrr::furrr_options(seed = NULL),
      # .progress = TRUE,
      db_path = db_path,
      perc_id = perc_id,
      perc_qcov_hsp = perc_qcov_hsp
    )
  } else {
    blast_res <- purrr::map_dfr(
      asvs,
      get_blast_results,
      num_thread = 1,
      db_path = db_path,
      perc_id = perc_id,
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
  if (!is.na(out_rds)) {
    readr::write_rds(
      x = blast_res,
      file = out_rds
    )
  }
  return(blast_res)
}

