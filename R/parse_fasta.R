#' @title Parsing fasta file to sequences only
#'
#' @description Extract sequences From FASTA file
#'
#' @param file_path Path to FASTA file.
#'
#' @returns vector containing the sequences.
#'
#' @examples
#' fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
#' seqs <- parse_fasta(file_path = fasta_path)
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
  final_seq_vector <- vector()
  current_seq_string <- character(0L)
  for (line in fasta_file) {
    if (startsWith(line, "#")) {
      final_seq_vector <- c(final_seq_vector, current_seq_string)
      current_seq_string <- character(0L)
      next
    }
    if (startsWith(line, ">")) {
      final_seq_vector <- c(final_seq_vector, current_seq_string)
      current_seq_string <- character(0L)
      next
    } else {
      current_seq_string <- paste0(current_seq_string, line, collapse = "")
    }
  }
  final_seq_vector <- c(final_seq_vector, current_seq_string)
  final_seq_vector <- final_seq_vector[final_seq_vector != ""]

  return(final_seq_vector)
}
