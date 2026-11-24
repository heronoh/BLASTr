#' Determine if a Command is Available and Install if Necessary
#'
#' Checks if a BLAST+ command-line tool (`blastn`, `makeblastdb`, ...) is available
#'   on the system.
#'   If not, it creates a conda environment and installs the required tool.
#'
#' @param cmd Character string specifying the command-line tool to check for.
#'   BLAST+ programs (`"blastn"`, `"tblastn"`, `"makeblastdb"`, ...). Default is `"blastn"`.
#' @param env_name Character string specifying the name of the conda
#'   environment to create or use.
#'   Default is `"blastr-blast-env"`.
#' @param verbose Character string specifying the verbosity level during
#'   environment creation.
#'   Options are `"silent"`, `"output"`, "cmd", or "full".
#'   Default is `"silent"`.
#' @param force Logical indicating whether to force the re-creation of the
#'   conda environment even if it exists.
#'   Default is `FALSE`.
#'
#' @returns Invisibly returns `TRUE` if the command is available
#'   or successfully installed.
#'
#' @examples
#' \dontrun{
#' # Check if 'blastn' command is available, and install it if not
#' check_cmd("blastn", env_name = "blastr-blast-env")
#'
#' # Force re-creation of the conda environment and re-install 'blastn'
#' check_cmd("blastn", env_name = "blastr-blast-env", force = TRUE)
#' }
#'
#' @keywords internal
#' @noRd
check_cmd <- function(
  cmd = "blastn",
  env_name = "blastr-blast-env",
  verbose = c("silent", "cmd", "output", "full"),
  force = FALSE
) {
  verbose <- rlang::arg_match(verbose)
  packages_to_install <- NULL
  # supported_archs <- names(conda_packages_list)
  # Check if operational system is in the included list of available systems
  # if (isTRUE(stringr::str_detect(condathis::get_sys_arch(), pattern = "^Windows"))) {
  # }

  if (stringr::str_detect(string = cmd, pattern = "^(t?)blast|makeblastdb")) {
    packages_to_install <- "bioconda::blast==2.16"
  } else {
    cli::cli_abort(
      message = c(
        `x` = "Unsupported command: only BLAST+ programs ('blastn', 'tblastn', 'makeblastdb', ...) are supported."
      ),
      class = "blastr_check_cmd_unsupported_cmd"
    )
  }

  if (isTRUE(force)) {
    condathis::create_env(
      packages = packages_to_install,
      channels = c("conda-forge", "bioconda"),
      env_name = env_name,
      verbose = verbose,
      overwrite = force
    )
  } else if (isFALSE(condathis::env_exists(env_name))) {
    condathis::create_env(
      packages = packages_to_install,
      channels = c("conda-forge", "bioconda"),
      env_name = env_name,
      verbose = verbose
    )
  }
  # else if (isFALSE(condathis::env_exists(env_name))) {
  #  condathis::create_env(
  #    env_name = env_name,
  #    verbose = verbose
  #  )
  # }
  invisible(TRUE)
}
