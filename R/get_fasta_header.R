#' @title Get header from SubjectID
#'
#' @description Retrieve complete sequence header from BLAST DB based on SubjectIDs
#'
#' @inheritParams run_blast
#'
#' @param id SubjectID from BLAST results or any NCBI Nucleotide identifier.
#'
#' @return Complete identifier for the SubjectID as is on the database.
#'
#' @examples
#' \dontrun{
#' BLASTr::get_fasta_header(id = "AP011979.1", db_path = "/data/databases/nt/nt")
#' }
#' @export
get_fasta_header <- function(id,
                             db_path,
                             env_name = "blast-env") {
  if (is.null(db_path)) {
    db_path <- getOption(
      "BLASTr.db_path",
      default = NULL
    )
  }
  if (is.null(db_path)) {
    cli::cli_abort(
      message = "No BLAST database provided.",
      class = "blastr_missing_blast_db"
    )
  }

  blastdbcmd_res <- condathis::run(
    "blastdbcmd",
    "-db", db_path,
    "-entry", id,
    "-outfmt", "%t",
    env_name = env_name
  )

  return(blastdbcmd_res$stdout)
}
