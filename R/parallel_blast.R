#' @title Run Parallelized BLAST
#'
#' @description Run parallel BLAST for a set of sequences
#'
#' @inheritParams get_blast_results
#'
#' @param asvs Character vector with sequences
#' @param out_file Complete path to output .csv file on an existing folder.
#' @param out_RDS Complete path to output RDS file on an existing folder.
#' @param total_cores Total available cores to run BLAST in parallel.
#'   Check your max with *future::availableCores()*.
#' @param blast_type BLAST+ executable to be used on search.
#'
#' @return A tibble with the BLAST tabular output.
#'
#' @examples
#' \dontrun{
#' blast_res <- BLASTr::parallel_blast(
#'   asvs = ASVs_test, # vector of sequences to be searched
#'   db_path = "/data/databases/nt/nt", # path to a formatted blast database
#'   out_file = NULL, # path to a .csv file to be created with results (on an existing folder)
#'   out_RDS = NULL, # path to a .RDS file to be created with results (on an existing folder)
#'   perc_id = 80, # minimum identity percentage cutoff
#'   perc_qcov_hsp = 80, # minimum percentage coverage of query sequence by subject sequence cutoff
#'   num_threads = 1, # number of threads/cores to run each blast on
#'   total_cores = 8, # number of total threads/cores to allocate all blast searches
#'   # maximum number of alignments/matches to retrieve results for each query sequence
#'   num_alignments = 3,
#'   blast_type = "blastn" # blast search engine to use
#' )
#' }
#' @export
parallel_blast <- function(
    asvs,
    db_path,
    out_file = NULL,
    out_RDS = NULL,
    num_threads = 1L,
    blast_type = "blastn",
    total_cores = 1L,
    perc_id = 80L,
    perc_qcov_hsp = 80L,
    num_alignments = 4L,
    verbose = "silent",
    env_name = "blastr-blast-env") {
  rlang::check_required(asvs)
  rlang::check_required(db_path)

  check_cmd(blast_type, env_name = env_name, verbose = verbose)

  parallel_set <- FALSE

  if (
    isTRUE(total_cores > 1L) &&
      isTRUE(mirai::status(.compute = "blastr-cpu")$connections < total_cores)
  ) {
    mirai::daemons(n = total_cores, .compute = "blastr-cpu")
    parallel_set <- TRUE
  }

  blast_res <- purrr::map(
    .x = asvs,
    .f = purrr::in_parallel(
      .f = \(x) {
        get_blast_results(
          asv = x,
          db_path = db_path,
          num_threads = num_threads,
          blast_type = blast_type,
          perc_id = perc_id,
          perc_qcov_hsp = perc_qcov_hsp,
          num_alignments = num_alignments,
          verbose = verbose,
          env_name = env_name
        )
      },
      num_threads = 1L,
      blast_type = blast_type,
      num_alignments = num_alignments,
      db_path = db_path,
      perc_id = perc_id,
      perc_qcov_hsp = perc_qcov_hsp,
      verbose = verbose,
      env_name = env_name,
      get_blast_results = get_blast_results,
      run_blast = run_blast
    )
  ) |>
    purrr::list_rbind()

  if (identical(class(blast_res), "data.frame")) {
    class(blast_res) <- c("tbl_df", "tbl", "data.frame")
  }

  if (isTRUE(!is.na(out_file) && !is.null(out_file))) {
    readr::write_csv(
      x = blast_res,
      file = out_file,
      append = FALSE
    )
  }

  if (isTRUE(!is.na(out_RDS) && !is.null(out_RDS))) {
    readr::write_rds(
      x = blast_res,
      file = out_RDS
    )
  }

  if (isTRUE(total_cores > 1L) && isTRUE(parallel_set)) {
    mirai::daemons(n = 0L, .compute = "blastr-cpu")
  }

  return(blast_res)
}
