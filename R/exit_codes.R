#' @title Exit codes and standard error of the BLAST runs
#'
#' @description
#' This function extracts the exit codes and standard error messages from the
#' results of a BLAST search performed using the [parallel_blast()] function.
#'
#' @param blast_res Results of [parallel_blast()].
#' @returns A tibble with one row per sequence: `query_seq`, `exit_code` of
#'   the `blastn` run of its batch (`NA` when cancelled) and `stderr`.
#'
#' @examples
#' # needs BLAST+ (see install_dependencies())
#' if (nzchar(Sys.which("blastn"))) {
#'   fasta <- system.file("extdata", "minimal_db_blast.fasta", package = "BLASTr")
#'   db <- file.path(tempdir(), "minimal_db_blast")
#'   make_blast_db(fasta, db_path = db)
#'   seqs <- readLines(system.file("extdata", "seqs_test.txt", package = "BLASTr"))
#'   res <- parallel_blast(seqs, db)
#'   exit_codes(res)
#' }
#'
#' @family search functions
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
  return(metadata$exit_codes) # nolint: return_linter
}
