#' Install Dependencies on Conda Environments
#' @export
install_dependencies <- function(verbose = "silent") {
  check_cmd("blastn", env_name = "blast-env", verbose = verbose)
  check_cmd("efetch", env_name = "entrez-env", verbose = verbose)
}
