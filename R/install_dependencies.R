#' @title Install BLAST+ in a conda environment
#'
#' @description
#' Installs BLAST+ 2.17.0 (`blastn`, `makeblastdb`, ...) from bioconda into
#' the conda environment `blastr-blast-env`, through the condathis package.
#' BLASTr runs BLAST+ from this environment; without it, from the system,
#' if it is BLAST+ 2.16.0 or newer (see [has_blast()]). When the environment
#' exists with another version of BLAST+ (created by an older BLASTr), it is
#' recreated with BLAST+ 2.17.0. Searches never install or update BLAST+: only this
#' function does. The online taxonomy
#' functions ([get_tax_by_taxID()], [get_tax_by_name()]) query the NCBI over
#' HTTP and need no external tool.
#'
#' @param verbose A character string specifying the verbosity level during environment creation.
#'   Options include `"silent"`, `"cmd"`, `"output"`, or `"full"`. Default is `"silent"`.
#' @param force Recreate the conda environment even if it exists. Default
#'   is `FALSE`.
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
  install_blast_env(env_name = "blastr-blast-env", verbose = verbose, force = force)
  invisible(TRUE)
}
