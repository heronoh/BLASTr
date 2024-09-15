# #' Install BLAST+ CLI through conda
# #'
# #' Install BLAST+ command line tools in a conda environment.
# #'   Currently it don't works on ARM 64-bits MacOS systems.
# #'
# #' @param conda_env_name Name of the conda environment where
# #'   BLAST+ will be installed. Default is "blastr-blast-env".
# #'
# #' @export
# install_blast_conda <- function(conda_env_name = "blastr-blast-env") {
#   if (isTRUE(Sys.which("blastn") != "")) {
#     cli::cli_inform(
#       "BLAST+ CLI is already installed. Skipping installation."
#     )
#     return()
#   }
#   if (
#     (Sys.info()["machine"] %in% "arm64") &&
#       (Sys.info()["sysname"] %in% "Darwin")
#   ) {
#     cli::cli_abort(
#       c(
#       "BLAST+ CLI for MacOS ARM 64-bits architecture should be installed with
#         Hombrew.",
#       i = "Use: {.code brew install blast}",
#       x = "conda recipe is not available for MacOS ARM 64-bits architecture,"
#       ),
#       class = "blastr_conda_package_not_available_platform"
#     )
#   }
#   if (!requireNamespace("reticulate", quietly = TRUE)) {
#     cli::cli_abort(
#       c(
#       x = "Package {.pkg reticulate} is required to install BLAST+ CLI through # conda.",
#       i = "Please install with {.code install.packages(\"reticulate\")}."
#     ),
#     class = "blastr_reticulate_not_installed"
#     )
#   }
#   reticulate::conda_install(
#     envname = "blastr-env",
#     packages = "blast",
#     channel = "bioconda"
#   )
# }
#
