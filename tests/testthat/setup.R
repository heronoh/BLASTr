withr::local_options(
  .new = list(
    warnPartialMatchArgs = TRUE,
    warnPartialMatchAttr = TRUE,
    warnPartialMatchDollar = TRUE
  ),
  .local_envir = testthat::teardown_env()
)

# The conda environment of BLASTr (with BLAST+) is the one already installed:
# keep its data directory when HOME and the other user directories are moved
# to temporary ones below.
conda_data_home <- fs::path_dir(fs::path_dir(fs::path_dir(fs::path_dir(conda_env_dir("x")))))

tmp_home_path <- withr::local_tempdir(
  pattern = "tmp-home",
  .local_envir = testthat::teardown_env()
)
tmp_data_path <- withr::local_tempdir(
  pattern = "tmp-data",
  .local_envir = testthat::teardown_env()
)
tmp_cache_path <- withr::local_tempdir(
  pattern = "tmp-cache",
  .local_envir = testthat::teardown_env()
)

tmp_wd_path <- fs::path(tmp_home_path, "wd")

if (isFALSE(fs::dir_exists(tmp_home_path))) {
  fs::dir_create(tmp_home_path)
}
if (isFALSE(fs::dir_exists(tmp_data_path))) {
  fs::dir_create(tmp_data_path)
}
if (isFALSE(fs::dir_exists(tmp_cache_path))) {
  fs::dir_create(tmp_cache_path)
}
if (isFALSE(fs::dir_exists(tmp_wd_path))) {
  fs::dir_create(tmp_wd_path)
}

withr::local_dir(
  new = tmp_wd_path,
  .local_envir = testthat::teardown_env()
)

withr::local_envvar(
  .new = list(
    `HOME` = tmp_home_path,
    `USERPROFILE` = tmp_home_path,
    `LOCALAPPDATA` = tmp_data_path,
    `APPDATA` = tmp_data_path,
    `R_USER_DATA_DIR` = conda_data_home,
    `XDG_DATA_HOME` = conda_data_home,
    `R_USER_CACHE_DIR` = tmp_cache_path,
    `XDG_CACHE_HOME` = tmp_cache_path
  ),
  .local_envir = testthat::teardown_env()
)

# `tmp_blast_db_path` is used in the actual tests.
# `make_blast_db()` writes the database next to the FASTA file,
# so the FASTA is copied to the temporary directory first.
tmp_blast_db_path <- fs::path(tmp_wd_path, "minimal_db_blast", ext = "fasta")
base::rm(tmp_home_path)
base::rm(tmp_data_path)
base::rm(tmp_cache_path)
base::rm(tmp_wd_path)
base::rm(conda_data_home)


# The tests that run BLAST+ use the conda environment of BLASTr as it is
# (they are skipped without it; see helper-blast.R): nothing is installed here.
if (isTRUE(has_blast())) {
  fs::file_copy(
    path = fs::path_package(
      "BLASTr",
      "extdata",
      "minimal_db_blast",
      ext = "fasta"
    ),
    new_path = tmp_blast_db_path,
    overwrite = TRUE
  )
  make_blast_db(
    fasta_path = tmp_blast_db_path,
    db_type = "nucl",
    verbose = "silent"
  )
  testthat::expect_true(
    fs::file_exists(fs::path(tmp_blast_db_path, ext = "nsq"))
  )
}
