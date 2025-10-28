testthat::test_that("parallel_blast works", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  ASVs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "asvs_test",
    ext = "txt"
  ))
  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast(
    asvs = ASVs_test,
    db_path = db_path,
    out_file = NA,
    out_RDS = NA,
    total_cores = 4,
    perc_id = 80,
    num_threads = 1,
    perc_qcov_hsp = 80,
    num_alignments = 4,
    blast_type = "blastn",
    verbose = "silent"
  )

  testthat::expect_s3_class(blast_res, "tbl_df")
})

testthat::test_that("run_blast works", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  ASVs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "asvs_test",
    ext = "txt"
  ))
  db_path <- tmp_blast_db_path

  blast_res <- run_blast(
    asv = ASVs_test[1],
    db_path = db_path,
    verbose = "silent"
  )

  testthat::expect_equal(blast_res$status, 0)
})

testthat::test_that("run_blast fails", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  ASVs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "asvs_test",
    ext = "txt"
  ))
  # `tmp_blast_db_path <- fs::file_temp("minimal_db_blast_")`
  db_path <- tmp_blast_db_path

  testthat::expect_error(
    blast_res <- run_blast(
      asv = ASVs_test[1],
      db_path = db_path,
      perc_id = "MISSING_VALUE_STRInG",
      verbose = "silent"
    ),
    class = "blastr_error_blast_run"
  )
  # testthat::expect_equal(blast_res$status, 0)
})
