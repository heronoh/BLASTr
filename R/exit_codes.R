#' Retrieve Exit Codes and Standard Error from BLASTr Results
#'
#' @description
#' This function extracts the exit codes and standard error messages from the
#' results of a BLAST search performed using the [parallel_blast()] function.
#'
#' @param blast_res A data frame containing the results of a BLAST search.
#' @return A data frame with the exit codes and standard error messages.
#' @export
exit_codes <- function(blast_res) {
  if (isFALSE("BLASTr_metadata" %in% names(attributes(blast_res)))) {
    cli::cli_abort(
      message = c(
        `x` = "The provided results do not contain {.field BLASTr_metadata}."
      ),
      class = "blastr_missing_metadata_error"
    )
  }
  metadata <- attributes(blast_res)$BLASTr_metadata

  if (is.null(metadata$exit_codes)) {
    cli::cli_abort(
      message = c(
        `x` = "No exit codes found in the provided BLAST results metadata."
      ),
      class = "blastr_missing_exit_codes_error"
    )
  }
  return(metadata$exit_codes)
}