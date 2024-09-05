#' @title Run BLASTn
#'
#' @description Retrieve BLAST results from a given ASV
#'
#' @inheritParams run_blast
#'
#' @return A `tibble` with the results of BLASTn for each sequence.
#'
#' @examples
#' blast_res <- BLASTr::get_blast_results(
#'   asv = "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTA
#'   CTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
#'   db_path = "/data/databases/nt/nt",
#'   perc_id = 80,
#'   num_thread = 1,
#'   perc_qcov_hsp = 80,
#'   num_alignments = 2
#' )
#'
#' @export
get_blast_results <- function(asv,
                              db_path,
                              num_threads,
                              blast_type = "blastn",
                              perc_id = 80,
                              perc_qcov_hsp = 80,
                              num_alignments = 4) {
  #  if (is.null(db_path)) {
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
  if (is.null(num_threads)) {
    num_threads <- getOption(
      "BLASTr.num_threads",
      default = 1
    )
  }

  blast_res <- run_blast(
    asv = asv,
    blast_type = blast_type,
    num_threads = num_threads,
    db_path = db_path,
    perc_id = perc_id,
    perc_qcov_hsp = perc_qcov_hsp,
    num_alignments = num_alignments
  )
  # test blast_res content ----
  if (blast_res$status != 0) {
    cli::cli_abort(
      message = "{.pkg BLASTr} has not run correctly."
    )
  }
  if (blast_res$stdout == "") {
    df_to_return <- tibble::tibble(`OTU` = asv)
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
        "staxid",
        show_col_types = FALSE
      ),
      trim_ws = TRUE,
      comment = "#"
      ) |>
    dplyr::mutate("staxid" = as.character(staxid)) # adicionado para testes

    # dplyr::mutate(dplyr::across(.cols = dplyr::ends_with("taxid"), # adicionado para testes
                         # ~ as.character(.)))

  blast_table$`subject header` <- purrr::map_chr(
    .x = blast_table$subject,
    .f = get_fasta_header,
    db_path = db_path
  )
  blast_table <- dplyr::relocate(
    blast_table,
    "subject header"
  )
  blast_table <- tibble::rowid_to_column(
    blast_table,
    var = "res"
  )
  # blast_table <- tidyr::pivot_wider(blast_table, id_cols = subject,  names_from = res, values_from = seq_len(ncol(blast_table)), names_glue = "{res}_{.value}")
  blast_table <- tidyr::pivot_wider(
    blast_table,
    names_from = "res",
    values_from = base::seq_len(
      base::ncol(blast_table)
    ),
    names_glue = "{res}_{.value}"
  )

  blast_table <- blast_table |>
    dplyr::mutate(`Sequence` = asv) |>
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
