#' Install Dependencies on Conda Environments
#' @export
install_dependencies <- function() {
  check_cmd("blastn", env_name = "blast-env")
  check_cmd("efetch", env_name = "entrez-env")
}
