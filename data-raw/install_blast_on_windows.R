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
        channels = c("conda-forge"),
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

# install_blast_windows()

check_cmd("blastn", env_name = "blastr-blast-win-env")
check_cmd("efetch", env_name = "blastr-blast-win-env")
check_cmd("makeblastdb", env_name = "blastr-blast-win-env")

# install_dependencies <- function(env_name = "blast-env") {
#   condathis::create_env("bioconda::entrez-direct==24.0", env_name = "entrez-env")
#
#  condathis::create_env("bioconda::blast==2.16", env_name = "blastblast-env")
# }

# Testing on Kappa
# ssh kappa
# rig add 4.5.1
# rig default 4.5.1
# R -q -s -e "message(R.version.string)"
#> R version 4.5.1 (2025-06-13 ucrt)
# dir Documents\projects
# cd Documents\projects
# R
# install.packages("pak")
# pak::pak("condathis")
# library(condathis)
# pak::pak("heronoh/BLASTr")
# pak::pak("luciorq/condathis")
# pak::pak("tictoc")

condathis::install_micromamba(force = TRUE)

# install_blast_windows()

library(BLASTr)
(BLASTr:::check_cmd("blastn", env_name = "blastr-blast-win-env", verbose = "full"))
(BLASTr:::check_cmd("efetch", env_name = "blastr-blast-win-env"))
(BLASTr:::check_cmd("makeblastdb", env_name = "blastr-blast-win-env"))


fasta_path <- fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")

# Path to database
db_path <- fs::path_temp("minimal_db_blast")

make_blast_db(
  fasta_path = fasta_path,
  db_path = db_path,
  db_type = "nucl",
  verbose = "full",
  env_name = "blastr-blast-win-env"
)

asvs <- c(
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGCGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "CTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAGAGTACTACAAGTGAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC",
  "GCCAAATTTGTGTTTTGTCCTTCGTTTTTAGTTAATTGTTACTGGCAAATGACTAACGACAAATGATAAATTACTAATAC",
  "AACATTGTATTTTGTCTTTGGGGCCTGGGCAGGTGCAGTAGGAACTTCACTTAGAATAATTATTCGTACTGAGCTTGGGCATCCAGGAAGACTTATCGGGGATGATCAAATCTATAATGTAATTGTTACAGCACATGCATTTGTGATAATTTTTTTTATAGTAATACCTATTATGATT",
  "ACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCCATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGGTGCCCCCGATATG",
  "TTAGCCATAAACATAAAAGTTCACATAACAAGAACTTTTGCCCGAGAACTACTAGCAACAGCTTAAAACTCAAAGGACTTGGCGGTGCTTTATATCCAC"
)


#> Parallel BLAST 8 processes with 3 threads each: 13.14 sec elapsed
tictoc::tic("Parallel BLAST 8 processes with 3 threads each")
res <- parallel_blast(
  asv = rep(asvs, 50),
  db_path = db_path,
  total_cores = 8L,
  num_threads = 3L,
  verbose = "silent",
  env_name = "blastr-blast-win-env"
)
tictoc::toc()

#> Parallel BLAST 8 processes with 2 threads each: 12.55 sec elapsed
tictoc::tic("Parallel BLAST 8 processes with 2 threads each")
res <- parallel_blast(
  asv = rep(asvs, 50),
  db_path = db_path,
  total_cores = 8L,
  num_threads = 2L,
  verbose = "silent",
  env_name = "blastr-blast-win-env"
)
tictoc::toc()

#> Parallel BLAST 8 processes with 1 threads each: 12.97 sec elapsed
tictoc::tic("Parallel BLAST 8 processes with 1 threads each")
res <- parallel_blast(
  asv = rep(asvs, 50),
  db_path = db_path,
  total_cores = 8L,
  num_threads = 1L,
  verbose = "silent",
  env_name = "blastr-blast-win-env"
)
tictoc::toc()

#> Parallel BLAST 2 processes with 2 threads each: 13.72 sec elapsed
tictoc::tic("Parallel BLAST 2 processes with 2 threads each")
res <- parallel_blast(
  asv = rep(asvs, 50),
  db_path = db_path,
  total_cores = 2L,
  num_threads = 2L,
  verbose = "silent",
  env_name = "blastr-blast-win-env"
)
tictoc::toc()

# Parallel BLAST 1 processes with 2 threads each: 16.18 sec elapsed
tictoc::tic("Parallel BLAST 1 processes with 2 threads each")
res <- parallel_blast(
  asv = rep(asvs, 50),
  db_path = db_path,
  total_cores = 1L,
  num_threads = 2L,
  verbose = "silent",
  env_name = "blastr-blast-win-env"
)
tictoc::toc()

tictoc::tic("Parallel BLAST 1 processes with 1 threads each")
res <- parallel_blast(
  asv = rep(asvs, 50),
  db_path = db_path,
  total_cores = 1L,
  num_threads = 1L,
  verbose = "silent",
  env_name = "blastr-blast-win-env"
)
tictoc::toc()

condathis::get_env_dir("blastr-blast-win-env") |>
  fs::dir_ls(recursive = TRUE)

parallel_get_tax(
  organisms_taxIDs = c("9606", "562", "3702", "9913", "10090", "10116"),
  total_cores = 1L,
  verbose = "full",
  env_name = "blastr-blast-win-env",
)

(BLASTr:::check_cmd("efetch", env_name = "blastr-blast-win-env", verbose = "full"))

get_tax_by_taxID(
  organisms_taxIDs = "9606",
  env_name = "blastr-blast-win-env",
  verbose = "full"
) |>
  print()

condathis::run(
  "efetch", "-help",
  env_name = "blastr-blast-win-env",
  verbose = "full"
)

condathis::get_env_dir("blastr-blast-win-env") |>
  fs::path("Library", "bin") |>
  fs::dir_ls() # recursive = TRUE)


condathis::run_bin("blastn", "-help", env_name = "blastr-blast-win-env", verbose = "full")

condathis::run("blastn", "-version", env_name = "blastr-blast-win-env", verbose = "full")
condathis::run_bin("blastn", "-version", env_name = "blastr-blast-win-env", verbose = "full")


condathis::run_bin("perl", "-help", env_name = "blastr-blast-win-env", verbose = "full")


condathis::run("where", "blastn", env_name = "blastr-blast-win-env", verbose = "full")

condathis::run("where", "R", env_name = "blastr-blast-win-env", verbose = "full")

condathis::run("where", "Rscript", env_name = "blastr-blast-win-env", verbose = "full")
