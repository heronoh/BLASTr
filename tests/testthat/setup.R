withr::local_options(
  .new = list(
    warnPartialMatchArgs = TRUE,
    warnPartialMatchAttr = TRUE,
    warnPartialMatchDollar = TRUE
  ),
  .local_envir = testthat::teardown_env()
)

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
    `R_USER_DATA_DIR` = tmp_data_path,
    `XDG_DATA_HOME` = tmp_data_path,
    `R_USER_CACHE_DIR` = tmp_cache_path,
    `XDG_CACHE_HOME` = tmp_cache_path
  ),
  .local_envir = testthat::teardown_env()
)

# `tmp_blast_db_path` is used in the actual tests.
tmp_blast_db_path <- fs::path(tmp_wd_path, "minimal_db_blast")
base::rm(tmp_home_path)
base::rm(tmp_data_path)
base::rm(tmp_cache_path)
base::rm(tmp_wd_path)


# This is the check used inside skip_if_cran
if (!(!interactive() && !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false"))))) {
  cli::cli_inform(c(`!` = "Setting up test environment for {.pkg BLASTr} ..."))

  testthat::expect_true(install_dependencies(verbose = "silent", force = FALSE))

  testthat::expect_true(check_cmd("blastn", env_name = "blastr-blast-env"))

  testthat::expect_true(check_cmd("tblastn", env_name = "blastr-blast-env"))

  testthat::expect_true(check_cmd("efetch", env_name = "blastr-entrez-env"))

  testthat::expect_true(check_cmd("makeblastdb", env_name = "blastr-blast-env"))

  # For local testing use:
  # + `tmp_blast_db_path <- fs::file_temp("minimal_db_blast_")`
  make_blast_db(
    fasta_path = fs::path_package(
      "BLASTr",
      "extdata",
      "minimal_db_blast",
      ext = "fasta"
    ),
    db_path = tmp_blast_db_path,
    db_type = "nucl",
    verbose = "silent"
  )

  testthat::expect_true(
    fs::file_exists(fs::path(tmp_blast_db_path, ext = "nsq"))
  )
}
