# ' Determine if Command is Available
#'
#' Check if Commands Are Available otherwise create conda environments for each tool
check_cmd <- function(cmd = "blastn", env_name = "blast-env", verbose = "silent", force = FALSE) {
  cmd_bin <- Sys.which(cmd)

  if (stringr::str_detect(cmd, "blast")) {
    packages_to_install <- "bioconda::blast==2.16"
  } else if (stringr::str_detect(cmd, "efetch")) {
    packages_to_install <- "bioconda::entrez-direct==22.4"
  }

  if(isTRUE(force)) {
    condathis::create_env(packages_to_install, env_name = env_name, verbose = verbose, overwrite = force)
  }

  if (isFALSE(condathis::env_exists(env_name)) && !nzchar(cmd_bin)) {
    condathis::create_env(packages_to_install, env_name = env_name, verbose = verbose)
  } else if (!condathis::env_exists(env_name) && nzchar(cmd_bin)) {
    condathis::create_env(env_name = env_name, verbose = verbose)
  }
  invisible(TRUE)
}
