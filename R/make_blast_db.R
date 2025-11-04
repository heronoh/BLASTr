#' @title Make BLAST Database
#'
#' @description Create a BLAST database from a FASTA file.
#'
#' @param fasta_path Path to the input FASTA file.
#' @param db_path Path to the output BLAST database.
#' @param db_type Type of database to create, either "nucl" or "prot".
#' @param taxid_map Optional path to a file mapping sequence IDs to
#'  taxonomy IDs.
#' @param parse_seqids Whether to parse sequence IDs.
#' @param verbose Should `[condathis::run()]` internal command be shown?
#' @param env_name The name of the conda environment with the parameter
#' (i.e. "blastr-blast-env").
#' @return The result of the `makeblastdb` command.
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
#' }
#'
#' @export
make_blast_db <- function(
  fasta_path,
  db_path,
  db_type = "nucl",
  taxid_map = NULL,
  parse_seqids = TRUE,
  verbose = "silent",
  env_name = "blastr-blast-env"
) {
  rlang::check_required(fasta_path)
  rlang::check_required(db_path)
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
