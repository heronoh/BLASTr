#' Check that a BLAST+ program is available
#'
#' Checks that a BLAST+ command-line tool (`blastn`, `makeblastdb`, ...) of
#'   version `blast_version_min` or newer is available: in the conda
#'   environment of BLASTr or, without it, on the system (see blast_bin.R).
#'   Nothing is installed: without BLAST+ it stops and says to run
#'   install_dependencies(). Tells once per session when the BLAST+ used is
#'   not the version installed by install_dependencies().
#'
#' @param cmd BLAST+ program (`"blastn"`, `"tblastn"`, `"makeblastdb"`, ...).
#' @param env_name Name of the conda environment of BLASTr.
#' @param verbose Ignored (kept for the callers).
#'
#' @returns `TRUE`, invisibly, or an error of class
#'   `blastr_error_blast_not_found`.
#'
#' @keywords internal
#' @noRd
check_cmd <- function(
  cmd = "blastn",
  env_name = "blastr-blast-env",
  verbose = c("silent", "cmd", "output", "full")
) {
  verbose <- rlang::arg_match(verbose)
  if (!stringr::str_detect(string = cmd, pattern = "^(t?)blast|makeblastdb")) {
    cli::cli_abort(
      message = c(
        `x` = "Unsupported command: only BLAST+ programs ('blastn', 'tblastn', 'makeblastdb', ...) are supported."
      ),
      class = "blastr_check_cmd_unsupported_cmd"
    )
  }
  b <- resolve_blast_bin(cmd, env_name)
  key <- paste(b$path, b$version)
  if (!key %in% blastr_state$notified) {
    if (identical(b$source, "system")) {
      cli::cli_inform(c(
        `i` = "Using the BLAST+ {b$version} of the system ({.path {b$path}}): the conda environment {.val {env_name}} has no BLAST+ {blast_version_min} or newer.",
        `i` = "{.code install_dependencies()} installs BLAST+ {blast_version_pinned} there."
      ))
    } else if (!identical(b$version, blast_version_pinned)) {
      cli::cli_inform(c(
        `i` = "BLAST+ {b$version} in the conda environment {.val {env_name}}; BLASTr is tested with {blast_version_pinned}.",
        `i` = "{.code install_dependencies()} updates it."
      ))
    }
    blastr_state$notified <- c(blastr_state$notified, key)
  }
  invisible(TRUE)
}
