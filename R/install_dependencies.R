#' Install Required Command-line Tool Dependencies
#'
#' Installs the necessary command-line tools ('blast' and 'efetch') into conda environments if they are not already available. This ensures that all dependencies required by the package are present and up-to-date.
#'
#' @param verbose A character string specifying the verbosity level during environment creation.
#'   Options include `"silent"`, `"cmd"`, `"output"`, or `"full"`. Default is `"silent"`.
#' @param force A logical value indicating whether to force the re-creation of the Conda environments even if they already exist. Default is `FALSE`.
#'
#' @returns Invisibly returns `TRUE` after attempting to install the dependencies.
#'
#' @examples
#' \dontrun{
#' # Install dependencies with default settings
#' install_dependencies()
#'
#' # Install dependencies with verbose output
#' install_dependencies(verbose = "output")
#'
#' # Force re-installation of dependencies
#' install_dependencies(force = TRUE)
#' }
#' @export
install_dependencies <- function(
  verbose = c("silent", "cmd", "output", "full"),
  force = FALSE
) {
  verbose <- rlang::arg_match(verbose)
  check_cmd(
    cmd = "blastn",
    env_name = "blastr-blast-env",
    verbose = verbose,
    force = force
  )
  check_cmd(
    cmd = "efetch",
    env_name = "blastr-entrez-env",
    verbose = verbose,
    force = force
  )
  invisible(TRUE)
}
