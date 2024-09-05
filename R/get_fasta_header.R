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
#' BLASTr::get_fasta_header(id = "AP011979.1", db_path = "/data/databases/nt/nt")
#'
#' @export

get_fasta_header <- function(id,
                             db_path) {
  if (is.null(db_path)) {
    db_path <- getOption(
      "BLASTr.db_path",
      default = NULL
    )
  }
  if (is.null(db_path)) {
    cli::cli_abort(
      message = "No BLAST database provided."
    )
  }
  blastdbcmd_bin <- check_bin("blastdbcmd")
  command <- "{blastdbcmd_bin} -db {db_path} -entry {id} -outfmt %t"
  result <- shell_exec(
    cmd = command
  )
  return(result$stdout)
}
