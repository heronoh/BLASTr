# Create empty env


install_blast_windows <- function(env_name = "blast-win-env") {
  condathis::create_env(
    env_name = env_name
  )

  fs::dir_ls(
    path = fs::path(condathis::get_install_dir(), "envs", env_name)
  )


  env_path <- fs::path(condathis::get_install_dir(), "envs", env_name)
  env_bin_path <- fs::path(env_path, "bin")
  if (isFALSE(fs::dir_exists(env_bin_path))) {
    fs::dir_create(env_bin_path)
  }
  # env_name <- "blastr-blast-win-env"

  # For blast+ 2.16
  # + <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-win64.exe>
  # + <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-win64.exe>

  # download.file()



  cli::cli_inform(
    c(
      `v` = "{.field BLAST} successfully installed.",
      `!` = "Using path: {.path {env_path}}."
    )
  )
}



check_bin("blast")

install_dependencies <- function(env_name = "blast-env") {
  condathis::create_env("bioconda::entrez-direct==22.4", env_name = "entrez-env")

  condathis::create_env("bioconda::blast==2.16", env_name = "blast-env")
}
