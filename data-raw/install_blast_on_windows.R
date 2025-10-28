install_blast_windows <- function(
  env_name = "blastr-blast-win-env",
  verbose = "full"
) {
  # Create env with dependencies

  withr::with_envvar(
    new = list(
      CONDA_OVERRIDE_WIN = "1"
    ),
    code = {
      condathis::create_env(
        c("perl", "curl", "m2-bash"),
        platform = "win-64",
        env_name = env_name,
        verbose = verbose
      )
    },
    action = "replace"
  )


  fs::dir_ls(
    path = fs::path(condathis::get_install_dir(), "envs", env_name)
  )


  env_path <- fs::path(condathis::get_install_dir(), "envs", env_name)
  env_bin_path <- fs::path(env_path, "Library", "bin")
  if (isFALSE(fs::dir_exists(env_bin_path))) {
    fs::dir_create(env_bin_path)
  }

  # env_name <- "blastr-blast-win-env"

  # For blast+ 2.16
  # + download https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-win64.tar.gz
  # + unpack ncbi-blast-2.16.0+-x64-win64.tar.gz
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-win64.tar.gz",
    destfile = fs::path(env_bin_path, "ncbi-blast-2.16.0+-x64-win64.tar.gz")
  )
  # Extract the .tar.gz file into ~~bin~~ ncbi dir by default
  utils::untar(
    tarfile = fs::path(env_bin_path, "ncbi-blast-2.16.0+-x64-win64.tar.gz"),
    exdir = env_bin_path
  )

  # copy all content of fs::path(env_bin_path, "ncbi-blast-2.16.0+", "bin") to fs::path(env_bin_path)
  fs::file_copy(
    path = fs::dir_ls(
      path = fs::path(env_bin_path, "ncbi-blast-2.16.0+", "bin")
    ),
    new_path = env_bin_path,
    overwrite = TRUE
  )
  fs::file_delete(
    fs::path(env_bin_path, "ncbi-blast-2.16.0+-x64-win64.tar.gz")
  )
  fs::dir_delete(
    fs::path(env_bin_path, "ncbi-blast-2.16.0+")
  )
  # fs::dir_ls(
  #   path = fs::path(env_bin_path)
  # )
  fs::file_exists(
    fs::path(env_bin_path, "blastn.exe")
  )

  condathis::run(
    "blastn.exe", "-help",
    env_name = env_name,
    verbose = verbose
  )

  # For Entrez Direct CLI on msys2
  # + Dependencies: `curl`, `perl`, `bash`
  # + curl -O https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh
  # + sh ./install-edirect.sh

  # <https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz>
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz",
    destfile = fs::path(env_bin_path, "edirect.tar.gz")
  )

  utils::untar(
    tarfile = fs::path(env_bin_path, "edirect.tar.gz"),
    exdir = env_bin_path
  )
  fs::dir_ls(env_bin_path)

  cli::cli_inform(
    c(
      `v` = "{.field BLAST} successfully installed.",
      `!` = "Using path: {.path {env_path}}."
    )
  )
}

check_cmd("blastn", env_name = "blastr-blast-win-env")
check_cmd("efetch", env_name = "blastr-blast-win-env")
check_cmd("makeblastdb", env_name = "blastr-blast-win-env")
# install_dependencies <- function(env_name = "blast-env") {
#   condathis::create_env("bioconda::entrez-direct==24.0", env_name = "entrez-env")
#
#  condathis::create_env("bioconda::blast==2.16", env_name = "blastblast-env")
# }
