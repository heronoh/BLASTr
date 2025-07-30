#' Determine if a Command is Available and Install if Necessary
#'
#' Checks if a specified command-line tool ('blastn' or 'efetch') is available
#'   on the system.
#'   If not, it creates a conda environment and installs the required tool.
#'
#' @param cmd Character string specifying the command-line tool to check for.
#'   Supported commands are `"blastn"` and `"efetch"`. Default is `"blastn"`.
#' @param env_name Character string specifying the name of the conda
#'   environment to create or use.
#'   Default is `"blastr-blast-env"`.
#' @param verbose Character string specifying the verbosity level during
#'   environment creation.
#'   Options are `"silent"`, `"verbose"`, etc. Default is `"silent"`.
#' @param force Logical indicating whether to force the re-creation of the
#'   conda environment even if it exists.
#'   Default is `FALSE`.
#'
#' @return Invisibly returns `TRUE` if the command is available
#'   or successfully installed.
#'
#' @examples
#' \dontrun{
#' # Check if 'blastn' command is available, and install it if not
#' check_cmd("blastn", env_name = "blastr-blast-env")
#'
#' # Check if 'efetch' command is available, and install it if not
#' check_cmd("efetch", env_name = "blastr-entrez-env")
#'
#' # Force re-creation of the conda environment and re-install 'blastn'
#' check_cmd("blastn", env_name = "blastr-blast-env", force = TRUE)
#' }
check_cmd <- function(
    cmd = "blastn",
    env_name = "blastr-blast-env",
    verbose = "silent",
    force = FALSE) {
  packages_to_install <- NULL
  if (stringr::str_detect(cmd, "blast")) {
    packages_to_install <- "bioconda::blast==2.16"
  } else if (stringr::str_detect(cmd, "efetch")) {
    packages_to_install <- "bioconda::entrez-direct==24.0"
  } else {
    cli::cli_abort(
      message = c(
        `x` = "Unsupported command: only 'blast' variants and 'efetch' are supported."
      ),
      class = "blastr_check_cmd_unsupported_cmd"
    )
  }

  if (isTRUE(force)) {
    condathis::create_env(
      packages = packages_to_install,
      env_name = env_name,
      verbose = verbose,
      overwrite = force
    )
  } else if (isFALSE(condathis::env_exists(env_name))) {
    condathis::create_env(
      packages = packages_to_install,
      env_name = env_name,
      verbose = verbose
    )
  } else if (isFALSE(condathis::env_exists(env_name))) {
    condathis::create_env(
      env_name = env_name,
      verbose = verbose
    )
  }
  invisible(TRUE)
}
