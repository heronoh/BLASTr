parallel_blast2 <- function(
  query_seqs,
  db_path,
  ...,
  # out_file = NULL,
  # out_RDS = NULL,
  total_cores = 1L,
  num_threads = 1L,
  blast_type = "blastn",
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  num_alignments = 4L,
  retry_times = 10L,
  mt_mode = c("2", "1", "0"),
  verbose = "silent",
  env_name = "blastr-blast-env"
) {
  rlang::check_required(query_seqs)
  rlang::check_required(db_path)
  rlang::check_dots_empty()
  mt_mode <- rlang::arg_match(mt_mode)
  .data <- rlang::.data

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

  par_res <- purrr::map(
    .x = query_seqs,
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
      num_alignments = num_alignments,
      num_threads = num_threads,
      blast_type = blast_type,
      perc_id = perc_id,
      perc_qcov_hsp = perc_qcov_hsp,
      mt_mode = mt_mode,
      verbose = verbose,
      env_name = env_name,
      blast_cmd = blast_cmd
    )
  )

  blast_res_df <- par_res |>
    purrr::map(
      .f = \(x) {
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
              "subject header",
              "Sequence"

            ),
            show_col_types = FALSE,
            trim_ws = TRUE,
            comment = "#"
          ) |>
          tibble::rowid_to_column(var = "res")
      }
    ) |>
    purrr::list_rbind() |>
    dplyr::mutate("staxid" = as.character(.data$staxid)) |>
    dplyr::relocate("subject header", .after = "res")

  blast_res_df |>
    tidyr::pivot_wider(
      # id_cols = !dplyr::all_of("Sequence"),
      names_from = "res",
      values_from = -dplyr::all_of("Sequence"),
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
    dplyr::select(-dplyr::ends_with(c("_res", "_query")))
}
