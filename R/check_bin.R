#' Check CLI Tool
#'
#' Helper function to check if a CLI tool is available.
#'   This function fails if tool is not available on PATH.
#'
#' @param cmd Name of the tool to check.
#'
#' @return Path to tool executable, if available.
check_bin <- function(cmd) {
   cmd_bin <- Sys.which(cmd)
  if (isTRUE(cmd_bin == "")) {
    cli::cli_abort(
      message = c(
        "x" = "{.pkg {cmd}} is not available on PATH.",
        "i" = "Please install it and try again."
      )
    )
  }
  return(cmd_bin)
}
