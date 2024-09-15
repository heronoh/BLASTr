#' @title Run BLAST on parallel
#'
#' @description Run parallel BLAST for set of sequences
#'
#' @inheritParams run_blast
#'
#' @param out_file Complete path to output .csv file on an existing folder.
#' @param out_RDS Complete path to output RDS file on an existing folder.
#' @param total_cores Total available cores to parallelize BLAST. Chech your max with *future::availableCores()*
#' @param blast_type BLAST+ executable to be used on search.
#'
#' @return A tibble with the BLAST tabular output.
#'
#' @examples
#' \dontrun{
#' blast_res <- BLASTr::parallel_blast(
#'   asvs = ASVs_test, # vector of sequences to be searched
#'   db_path = "/data/databases/nt/nt", # path to a formated blast database
#'   out_file = NULL, # path to a .csv file to be created with results (on an existing folder)
#'   out_RDS = NULL, # path to a .RDS file to be created with results (on an existing folder)
#'   perc_id = 80, # minimum identity percentual cutoff
#'   perc_qcov_hsp = 80, # minimum percentual coverage of query sequence by subject sequence cutoff
#'   num_threads = 1, # number of threads/cores to run each blast on
#'   total_cores = 8, # number of tota threads/cores to alocate all blast searches
#'   num_alignments = 3, # maximum number of alignments/matches to retrieve results for each query sequence
#'   blast_type = "blastn" # blast search engine to use
#' )
#' }
#' @export

parallel_blast <- function(asvs,
                           db_path,
                           out_file,
                           out_RDS,
                           num_threads,
                           blast_type,
                           total_cores,
                           perc_id,
                           perc_qcov_hsp,
                           num_alignments) {
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
      workers = total_cores
    )

    blast_res <- furrr::future_map_dfr(
      .x = asvs,
      .f = get_blast_results,
      .options = furrr::furrr_options(seed = NULL),
      num_threads = 1,
      blast_type = blast_type,
      num_alignments = num_alignments,
      db_path = db_path,
      perc_id = perc_id,
      perc_qcov_hsp = perc_qcov_hsp
    )
  } else {
    blast_res <- purrr::map_dfr(
      .x = asvs,
      .f = get_blast_results,
      num_threads = 1,
      blast_type = blast_type,
      num_alignments = num_alignments,
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
  if (!is.na(out_RDS)) {
    readr::write_rds(
      x = blast_res,
      file = out_RDS
    )
  }
  return(blast_res)
}
