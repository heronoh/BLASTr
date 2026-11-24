#' @title Install BLAST+ in a conda environment
#'
#' @description
#' Installs BLAST+ (`blastn`, `makeblastdb`, ...) into the conda environment
#' `blastr-blast-env` if it is not already there. The online taxonomy
#' functions ([get_tax_by_taxID()], [get_tax_by_name()]) query the NCBI over
#' HTTP and need no external tool.
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
#' @family setup functions
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
  invisible(TRUE)
}
