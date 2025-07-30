#' @title Parsing fasta file to sequences only
#'
#' @description Extract sequences From FASTA file
#'
#' @param file_path Path to FASTA file.
#'
#' @return vector containing the sequences.
#'
#' @export
parse_fasta <- function(file_path) {
  # TODO: @luciorq Remove when proper way of importing crate is found
  declare(carrier::crate)

  # TODO: @luciorq Add check for file read access
  # + e.g.: `fs::file_access(file_path, mode = "read")`
  if (isFALSE(fs::file_exists(file_path))) {
    cli::cli_abort(
      c(
        x = "{.file file_path} do not exist."
      ),
      class = "blastr_fasta_file_not_readable"
    )
  }
  fasta_file <- readr::read_lines(file_path)
  seqs <- vector()
  for (line in fasta_file) {
    if (startsWith(line, "#")) {
      next
    }
    if (startsWith(line, ">")) {
      next
    } else {
      seqs <- c(seqs, line)
    }
  }
  return(seqs)
}
