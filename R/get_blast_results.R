#          function to get fasta names from db based on subjectIDs         ###

#' @title Run BLASTn
#'
#' @description Retrieve BLAST results from a given ASV
#'
#' @inheritParams run_blast
#'
#' @return A `tibble` with the results of BLASTn
#'
#' @export
<<<<<<< HEAD
get_blast_results <- function(
  asv,
  num_thread,
  total_cores,
  db_path,
  # blast_type,
  perc_ID,
  perc_qcov_hsp,
  num_alignments
) {
=======
get_blast_results <- function(asv,
                              num_thread,
                              total_cores,
                              db_path,
                              perc_ID,
                              perc_qcov_hsp) {
>>>>>>> be0da30911512133435167e5332b24284500284e
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
  if (is.null(num_thread)) {
    num_thread <- getOption(
      "BLASTr.num_thread",
      default = 1
    )
  }

  blast_res <- run_blast(
    asv = asv,
    # blast_type = blast_type,
    num_thread = num_thread,
    db_path = db_path,
    perc_ID = perc_ID,
    perc_qcov_hsp = perc_qcov_hsp,
    num_alignments = num_alignments
  )
  # test blast_res content ----
  if (blast_res$status != 0) {
    cli::cli_abort(
      message = "{.pkg BLASTr} has not run correctly."
    )
  }
  # `%>%` <- dplyr::`%>%`

  if (blast_res$stdout == "") {
    # rlang::inform(glue::glue("Sequence {asv} not found."))
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
        "qcovhsp"
      ),
      trim_ws = TRUE,
      comment = "#"
    )

  blast_table$`subject header` <- purrr::map_chr(
    .x = blast_table$subject,
    .f = get_fasta_header,
    db_path = db_path
  )
  blast_table <- dplyr::relocate(
    blast_table,
    `subject header`
  )
  blast_table <- tibble::rowid_to_column(
    blast_table,
    var = "res"
  )
  # blast_table <- tidyr::pivot_wider(blast_table, id_cols = subject,  names_from = res, values_from = seq_len(ncol(blast_table)), names_glue = "{res}_{.value}")
  blast_table <- tidyr::pivot_wider(
    blast_table,
    names_from = res,
    values_from = base::seq_len(
      base::ncol(blast_table)
    ),
    names_glue = "{res}_{.value}"
  )
  # TODO: @heron limpar colunas redundantes. Não está ativado pois estava dando pau
  # |>
  #     select(-c("1_res",
  #               "2_res",
  #               "3_res"))

  blast_table <- blast_table |>
    dplyr::mutate(`OTU` = asv) |>
    dplyr::relocate(tidyr::starts_with("6_")) |>
    dplyr::relocate(tidyr::starts_with("5_")) |>
    dplyr::relocate(tidyr::starts_with("4_")) |>
    dplyr::relocate(tidyr::starts_with("3_")) |>
    dplyr::relocate(tidyr::starts_with("2_")) |>
    dplyr::relocate(tidyr::starts_with("1_")) |>
    dplyr::relocate(`OTU`) |>
    dplyr::select(-tidyr::ends_with(c("_res","_query")))

  return(blast_table)
}
