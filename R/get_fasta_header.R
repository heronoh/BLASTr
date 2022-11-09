#'  Retrieve FASTA names from BLAST DB based on subjectIDs
#'
#' @param id SubjectID from BLAST results
#' @param db_path Path to BLAST DB
#'
#' @export
get_fasta_header <- function(
  id,
  db_path
) {
  # id <- "JQ365494.1"
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
