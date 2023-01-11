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
  fasta_file <- readr::read_lines(file_path)
  if (!fs::file_exists(file_path)) {
    cli::cli_abort(
      c(
        "{.file file_path} do not exist."
      )
    )
  }
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
