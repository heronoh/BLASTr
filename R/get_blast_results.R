#' @title Run BLASTn
#'
#' @description Retrieve BLAST results from a given ASV
#'
#' @inheritParams run_blast
#'
#' @returns A `tibble` with the results of BLASTn for each sequence.
#'
#' @examples
#' \dontrun{
#' dna_fasta_path <- fs::path_package(
#'   "BLASTr", "extdata", "minimal_db_blast",
#'   ext = "fasta"
#' )
#' temp_db_path <- fs::path_temp("minimal_db_blast")
#' make_blast_db(fasta_path = dna_fasta_path, db_path = temp_db_path)
#' asvs_string <- "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAG
#' AGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"
#'
#' get_blast_results(asv = asvs_string, db_path = temp_db_path)
#' }
#' @export
get_blast_results <- function(
  asv,
  db_path,
  num_threads = 1L,
  blast_type = "blastn",
  perc_id = 80L,
  perc_qcov_hsp = 80L,
  num_alignments = 4L,
  # task = task,
  # gapopen = 5,
  # gapextend = 2,
  verbose = c("silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env"
) {
  .data <- rlang::.data

  blast_res <- run_blast(
    asv = asv,
    blast_type = blast_type,
    num_threads = num_threads,
    db_path = db_path,
    perc_id = perc_id,
    perc_qcov_hsp = perc_qcov_hsp,
    num_alignments = num_alignments,
    verbose = verbose,
    env_name = env_name
  )

  if (isTRUE(blast_res$status != 0)) {
    cli::cli_abort(
      message = "{.pkg BLASTr} has not run correctly.",
      class = "blastr_run_blast_error"
    )
  }
  if (blast_res$stdout == "") {
    df_to_return <- tibble::tibble(`Sequence` = asv)
    return(df_to_return)
  }

  blast_table <- blast_res$stdout |>
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
        "staxid"
      ),
      show_col_types = FALSE,
      trim_ws = TRUE,
      comment = "#"
    ) |>
    dplyr::mutate("staxid" = as.character(.data$staxid))

  blast_table$`subject header` <- purrr::map_chr(
    .x = blast_table$subject,
    .f = get_fasta_header,
    db_path = db_path,
    verbose = verbose
  )

  blast_table <- dplyr::relocate(
    blast_table,
    "subject header"
  ) |>
    dplyr::mutate(
      `subject header` = stringr::str_trim(
        .data[["subject header"]],
        side = "both"
      )
    )

  blast_table <- tibble::rowid_to_column(
    blast_table,
    var = "res"
  )

  blast_table <- tidyr::pivot_wider(
    blast_table,
    names_from = "res",
    values_from = base::seq_len(
      base::ncol(blast_table)
    ),
    names_glue = "{res}_{.value}"
  )

  blast_table <- blast_table |>
    dplyr::mutate(`Sequence` = stringr::str_replace_all(asv, "\\s", "")) |>
    dplyr::relocate(tidyr::starts_with("6_")) |>
    dplyr::relocate(tidyr::starts_with("5_")) |>
    dplyr::relocate(tidyr::starts_with("4_")) |>
    dplyr::relocate(tidyr::starts_with("3_")) |>
    dplyr::relocate(tidyr::starts_with("2_")) |>
    dplyr::relocate(tidyr::starts_with("1_")) |>
    dplyr::relocate("Sequence") |>
    dplyr::select(-tidyr::ends_with(c("_res", "_query")))

  return(blast_table)
}
