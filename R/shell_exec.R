#' Run CLI Tools
#' @keywords internal
shell_exec <- function(cmd, .envir = parent.frame()) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    cli::cli_abort(
      message = "Package {.pkg processx} package is not installed."
    )
  }
  cmd_res <- processx::run(
    command = "bash",
    args = c(
      "-c",
      glue::glue(cmd, .envir = .envir)
    ),
    echo_cmd = FALSE
  )
  return(cmd_res)
}
