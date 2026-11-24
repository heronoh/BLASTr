#' @title Make BLAST Database
#'
#' @description Create a BLAST database from a FASTA file.
#'
#' @param fasta_path Path to the input FASTA file.
#' @param db_path Path (prefix) of the output BLAST database.
#'   Defaults to `fasta_path`, i.e. the database files are written next to the
#'   FASTA file (`<fasta_path>.nsq`, `<fasta_path>.nin`, ...).
#'   The function stops with an error if the directory of `db_path` does not
#'   exist or cannot be written to.
#' @param db_type Type of database to create, either "nucl" or "prot".
#' @param taxid_map Optional path to a file mapping sequence IDs to
#'   taxonomy IDs (two columns, separated by a space or a tab), passed to
#'   `makeblastdb -taxid_map`.
#' @param parse_seqids Whether to parse sequence IDs.
#' @param verbose Should `[condathis::run()]` internal command be shown?
#'   With `"silent"` (default) nothing is printed.
#' @param env_name The name of the conda environment with the parameter
#' (i.e. "blastr-blast-env").
#'
#' @returns The result of the `makeblastdb` command.
#'
#' @examples
#' \dontrun{
#' make_blast_db(
#'   fasta_path = fs::path_package(
#'     "BLASTr", "extdata", "minimal_db_blast",
#'     ext = "fasta"
#'   ),
#'   db_path = fs::path_temp("minimal_db_blast"),
#'   db_type = "nucl"
#' )
#'
#' # Without `db_path`, the database is written next to the FASTA file:
#' fasta_copy <- fs::path_temp("my_db", ext = "fasta")
#' fs::file_copy(
#'   fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta"),
#'   fasta_copy
#' )
#' make_blast_db(fasta_path = fasta_copy)
#' }
#'
#' @export
make_blast_db <- function(
  fasta_path,
  db_path = NULL,
  db_type = "nucl",
  taxid_map = NULL,
  parse_seqids = TRUE,
  verbose = c("silent", "cmd", "output", "full"),
  env_name = "blastr-blast-env"
) {
  rlang::check_required(fasta_path)
  db_path_missing <- rlang::is_null(db_path)
  if (isTRUE(db_path_missing)) {
    db_path <- fasta_path
  }
  check_db_dir_writable(
    db_path = db_path,
    fasta_path = fasta_path,
    db_path_missing = db_path_missing
  )
  verbose <- rlang::arg_match(verbose)
  check_cmd(cmd = "makeblastdb", env_name = env_name, verbose = verbose)
  args <- c(
    "-in",
    fasta_path,
    "-dbtype",
    db_type,
    "-out",
    db_path
  )

  if (isTRUE(parse_seqids)) {
    args <- c(args, "-parse_seqids")
  }

  if (!is.null(taxid_map)) {
    if (!isTRUE(fs::file_exists(taxid_map))) {
      cli::cli_abort(
        c(`x` = "{.arg taxid_map} does not exist: {.file {taxid_map}}"),
        class = "blastr_error_taxid_map_not_found"
      )
    }
    args <- c(args, "-taxid_map", taxid_map)
  }

  # NOTE: attemp to remove phone home in blast
  withr::local_envvar(
    .new = list(
      BLAST_USAGE_REPORT = "false"
    ),
    action = "replace"
  )
  blast_db_res <- condathis::run(
    cmd = "makeblastdb",
    args,
    env_name = env_name,
    verbose = verbose,
    error = "continue"
  )

  if (!identical(verbose, "silent")) {
    cli::cli_inform(c(`i` = "Command used: {.code makeblastdb {paste(args, collapse = ' ')}}"))
  }

  if (isTRUE(blast_db_res$status != 0L)) {
    error_msg_list <- stringr::str_extract_all(blast_db_res$stderr, "Error:.*")

    if (isTRUE(length(unlist(error_msg_list)) > 0)) {
      error_msg_vector <- unlist(error_msg_list)
      names(error_msg_vector) <- c(rep("x", times = length(error_msg_vector)))
    } else {
      error_msg_vector <- blast_db_res$stderr
    }

    cli::cli_abort(
      message = c(
        `!` = "Status code: {.val {blast_db_res$status}}",
        error_msg_vector
      ),
      class = "blastr_error_make_blast_db"
    )
  }

  return(invisible(blast_db_res))
}

#' Check that the BLAST database directory can be written to
#'
#' @param db_path Path (prefix) of the output BLAST database.
#' @param fasta_path Path to the input FASTA file.
#' @param db_path_missing Whether `db_path` was derived from `fasta_path`.
#' @param call Calling environment shown in the error.
#'
#' @returns `TRUE` invisibly, or an error of class
#'   `blastr_error_db_dir_not_writable`.
#'
#' @keywords internal
#' @noRd
check_db_dir_writable <- function(
  db_path,
  fasta_path,
  db_path_missing,
  call = rlang::caller_env()
) {
  db_dir <- fs::path_dir(fs::path_abs(db_path))
  dir_ok <- isTRUE(fs::dir_exists(db_dir)) &&
    isTRUE(unname(fs::file_access(db_dir, mode = "write")))
  if (isTRUE(dir_ok)) {
    return(invisible(TRUE))
  }
  reason <- if (isTRUE(fs::dir_exists(db_dir))) {
    "cannot be written to"
  } else {
    "does not exist"
  }
  if (isTRUE(db_path_missing)) {
    msg <- c(
      `x` = "The directory of {.arg fasta_path} {reason}: {.path {db_dir}}",
      `i` = "By default the BLAST database is written next to the FASTA file.",
      `>` = "Use {.arg db_path} to set where the formatted database will be written, e.g. {.code db_path = \"/a/writable/dir/{fs::path_file(fasta_path)}\"}."
    )
  } else {
    msg <- c(
      `x` = "The directory of {.arg db_path} {reason}: {.path {db_dir}}",
      `>` = "Choose a {.arg db_path} inside an existing directory with write permission."
    )
  }
  cli::cli_abort(
    msg,
    class = "blastr_error_db_dir_not_writable",
    call = call
  )
}
