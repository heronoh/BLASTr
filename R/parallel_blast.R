#' @title Run Parallelized BLAST
#'
#' @description Run parallel BLAST for a set of sequences
#'
#' @param query_seqs Character vector with sequences.
#' @param db_path Path to the formatted BLAST database.
#' @param total_cores Number of parallel BLAST processes to run.
#' @param num_threads Number of threads/cores to run each BLAST process on.
#' @param blast_type BLAST+ executable to be used on search.
#' @param perc_id Lowest identity percentage cutoff.
#' @param perc_qcov_hsp Lowest query coverage per HSP percentage cutoff.
#' @param num_alignments Number of alignments to retrieve results for each query
#'   sequence. Defaults to `4`.
#' @param retry_times Number of times to retry failed BLAST jobs.
#'   Defaults to `3` attempts.
#' @param mt_mode Multithreading mode to be used by BLAST+.
#'  One of: `c("2", "1", "0")`. See BLAST+ manual for details.
#' @param verbose Strategy used for showing outputting internal commands.
#'   Defaults to "progress".
#' @param env_name The name of the conda environment used to run
#'  command-line tools. Defaults to `"blastr-blast-env"`.
#'
#' @param asvs Deprecated. Same as `query_seqs`. Use `query_seqs` instead.
#' @param out_file Deprecated. Path to output `.csv` file on an existing directory.
#' @param out_RDS Deprecated. Path to output `RDS` file on an existing directory.
#'
#' @inheritParams rlang::args_dots_empty
#'
#' @return A tibble with the BLAST tabular output.
#'
#' @examples
#' \dontrun{
#' dna_fasta_path <- fs::path_package(
#'   "BLASTr", "extdata", "minimal_db_blast",
#'   ext = "fasta"
#' )
#' temp_db_path <- fs::path_temp("minimal_db_blast")
#' make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
#' asvs_string <- c(
#'   "CTAGCCATAAACTTAAATGAAGCTATACTAA",
#'   "ACTCGTTCGCCAGAGTACTACAAGCGAAAG"
#' )
#' blast_res <- parallel_blast(
#'   query_seqs = asvs_string,
#'   db_path = temp_db_path
#' )
#' blast_res
#' }
#'
#' @importFrom lifecycle deprecated
#'
#' @export
parallel_blast <- function(
  query_seqs,
  db_path,
  ...,
  total_cores = 1L,
  num_threads = 1L,
  blast_type = "blastn",
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  num_alignments = 4L,
  retry_times = 3L,
  mt_mode = c("2", "1", "0"),
  verbose = c("progress", "silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env",
  asvs = deprecated(),
  out_file = deprecated(),
  out_RDS = deprecated()
) {

  rlang::check_dots_empty()

  if (lifecycle::is_present(asvs)) {
    lifecycle::deprecate_warn(
      "0.1.7",
      "parallel_blast(asvs)",
      "parallel_blast(query_seqs)"
    )
    if (lifecycle::is_present(query_seqs)) {
      cli::cli_warn(
        "Both {.arg asvs} and {.arg query_seqs} were provided. Using {.arg query_seqs}."
      )
      asvs <- query_seqs
    }
    query_seqs <- asvs
  }
  if (lifecycle::is_present(out_file)) {
    lifecycle::deprecate_soft(
      "0.1.7",
      "parallel_blast(out_file)"
    )
  }
  if (lifecycle::is_present(out_RDS)) {
    lifecycle::deprecate_soft(
      "0.1.7",
      "parallel_blast(out_RDS)"
    )
  }

  rlang::check_required(query_seqs)
  rlang::check_required(db_path)
  mt_mode <- rlang::arg_match(mt_mode)

  # check for verbose
  verbose <- rlang::arg_match(verbose)
  if (
    rlang::is_interactive() &&
      isTRUE(verbose %in% c("progress", "full", "output"))
  ) {
    progress_var <- TRUE
  } else {
    progress_var <- FALSE
  }
  if (identical(verbose, "progress")) {
    verbose <- "silent"
  }

  # Check for query_seqs
  if (
    rlang::is_null(query_seqs) ||
      isTRUE(length(query_seqs) == 0L) ||
      rlang::is_na(query_seqs)
  ) {
    cli::cli_abort(
      c(
        `x` = "No query sequences were provided to {.fun parallel_blast}."
      ),
      class = "blastr_no_query_seqs_error"
    )
  }

  # Check for db_path
  if (
    rlang::is_null(db_path) ||
      isTRUE(length(db_path) == 0L) ||
      rlang::is_na(db_path)
  ) {
    cli::cli_abort(
      c(
        `x` = "{.pkg BLASTr}: No database path was provided to {.fun parallel_blast}."
      ),
      class = "blastr_no_db_path_error"
    )
  }
  .data <- rlang::.data
  .env <- rlang::.env

  check_cmd(blast_type, env_name = env_name, verbose = verbose)

  if (
    isTRUE(length(query_seqs) > 1L) &&
      isTRUE(total_cores > 1L) &&
      isTRUE(mirai::status()$connections < total_cores)
  ) {
    withr::defer(
      expr = {
        mirai::daemons(n = 0L)
      }
    )
    mirai::daemons(n = total_cores)
  }

  seqs_to_run <- base::unique(query_seqs)
  retry_count <- 0L
  query_seqs_final <- character(0L)
  par_res_final <- list()

  if (isTRUE(verbose %in% c("progress", "output", "full")) || !isFALSE(progress_var)) {
    cli::cli_inform(
      c(
        `i` = "Running {.fun parallel_blast} for {length(seqs_to_run)} unique query sequences.",
        `*` = "Using {total_cores} total processes with {num_threads} threads each."
      )
    )
  }
  while (
    isTRUE(retry_count <= retry_times) &&
      length(seqs_to_run) > 0L
  ) {
    par_res <- purrr::map(
      .x = seqs_to_run,
      .f = purrr::in_parallel(
        .f = \(x) {
          blast_cmd(
            query_str = x,
            db_path = db_path,
            num_alignments = num_alignments,
            num_threads = num_threads,
            blast_type = blast_type,
            perc_id = perc_id,
            perc_qcov_hsp = perc_qcov_hsp,
            mt_mode = mt_mode,
            verbose = verbose,
            env_name = env_name
          )
        },
        db_path = db_path,
        num_alignments = as.character(num_alignments),
        num_threads = num_threads,
        blast_type = blast_type,
        perc_id = perc_id,
        perc_qcov_hsp = perc_qcov_hsp,
        mt_mode = mt_mode,
        verbose = verbose,
        env_name = env_name,
        blast_cmd = blast_cmd
      ),
      .progress = progress_var
    )

    error_seqs <- par_res |>
      purrr::map2(
        .y = seqs_to_run,
        .f = \(x, y) {
          if (x[["status"]] != 0L) {
            return(y)
          } else {
            return(character(0L))
          }
        }
      ) |>
      unlist()

    retry_count <- retry_count + 1L
    par_res_final <- c(par_res_final, par_res)
    query_seqs_final <- c(query_seqs_final, seqs_to_run)
    seqs_to_run <- error_seqs
    if (
      isTRUE(length(seqs_to_run) > 0L) &&
        isTRUE(retry_count <= retry_times) &&
        (isTRUE(verbose %in% c("progress", "output", "full")) || !isFALSE(progress_var))
    ) {
      cli::cli_inform(
        c(
          `i` = "Retrying {.fun parallel_blast} for {length(seqs_to_run)} failed query sequences.",
          `*` = "Retry attempt {retry_count} of {retry_times}."
        )
      )
    }
  }

  blast_res_df <- par_res_final |>
    purrr::map2(
      .y = query_seqs_final,
      .f = \(x, y) {
        if (identical(x[["stdout"]], "")) {
          return(
            tibble::tibble(
              `res` = 1L,
              `query` = NA_character_,
              `subject` = NA_character_,
              # `indentity` = NA_real_,
              # `length` = NA_integer_,
              `staxid` = "N/A",
              `subject header` = NA_character_,
              `Sequence` = stringr::str_replace_all(y, "\\s", ""),
              `exit_code` = x[["status"]],
              `stderr` = x[["stderr"]]
            )
          )
        }
        x[["stdout"]] |>
          readr::read_delim(
            delim = "\t",
            col_names = c(
              "query",
              "subject",
              "indentity",
              "length",
              "mismatches",
              "gaps",
              "query start",
              "query end",
              "subject start",
              "subject end",
              "e-value",
              "bitscore",
              "qcovhsp",
              "staxid",
              "subject header"
            ),
            show_col_types = FALSE,
            trim_ws = TRUE,
            comment = "#"
          ) |>
          dplyr::mutate(
            `Sequence` = stringr::str_replace_all(.env$y, "\\s", ""),
            `exit_code` = x[["status"]],
            `stderr` = x[["stderr"]]
          ) |>
          tibble::rowid_to_column(var = "res")
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::mutate("staxid" = as.character(.data$staxid)) |>
    dplyr::relocate("subject header", .after = "res") |>
    dplyr::distinct()

  if (identical(class(blast_res_df), "data.frame")) {
    class(blast_res_df) <- c("tbl_df", "tbl", "data.frame")
  }

  blast_res <- blast_res_df |>
    tidyr::pivot_wider(
      names_from = "res",
      values_from = -dplyr::all_of(c("Sequence", "exit_code", "stderr")),
      names_glue = "{res}_{.value}",
      values_fn = list
    ) |>
    tidyr::unnest(cols = dplyr::everything()) |>
    dplyr::relocate(dplyr::starts_with("6_")) |>
    dplyr::relocate(dplyr::starts_with("5_")) |>
    dplyr::relocate(dplyr::starts_with("4_")) |>
    dplyr::relocate(dplyr::starts_with("3_")) |>
    dplyr::relocate(dplyr::starts_with("2_")) |>
    dplyr::relocate(dplyr::starts_with("1_")) |>
    dplyr::relocate("Sequence") |>
    dplyr::select(
      -dplyr::ends_with(c("_res", "_query")),
      -dplyr::starts_with("NA_")
    )

  attributes(blast_res)$BLASTr_metadata <- list(
    db_path = db_path,
    num_alignments = num_alignments,
    perc_id = perc_id,
    perc_qcov_hsp = perc_qcov_hsp,
    blast_type = blast_type,
    exit_codes = tibble::tibble(
      query_seq = blast_res_df$Sequence,
      exit_code = blast_res_df$exit_code,
      stderr = blast_res_df$stderr
    )
  )

  blast_res <- blast_res |>
    dplyr::select(-dplyr::any_of(c("exit_code", "stderr")))

  if (
    !rlang::is_na(out_file) &&
      !rlang::is_null(out_file) &&
      !rlang::is_missing(out_file)
  ) {
    readr::write_csv(
      x = blast_res,
      file = out_file,
      append = FALSE
    )
  }

  if (
    !rlang::is_na(out_RDS) &&
      !rlang::is_null(out_RDS) &&
      !rlang::is_missing(out_RDS)
  ) {
    readr::write_rds(
      x = blast_res,
      file = out_RDS
    )
  }

  return(blast_res)
}
