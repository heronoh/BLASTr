testthat::test_that("run_blast works", {
  skip_if_no_blast()

  seqs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "seqs_test",
    ext = "txt"
  ))
  db_path <- tmp_blast_db_path

  blast_res <- run_blast(
    query_seqs = seqs_test[1],
    db_path = db_path,
    verbose = "silent"
  )

  testthat::expect_equal(blast_res$status, 0)
})

testthat::test_that("run_blast fails", {
  skip_if_no_blast()

  seqs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "seqs_test",
    ext = "txt"
  ))
  db_path <- tmp_blast_db_path

  testthat::expect_error(
    run_blast(
      query_seqs = seqs_test[1],
      db_path = db_path,
      perc_id = "MISSING_VALUE_STRInG",
      verbose = "silent"
    ),
    class = "blastr_error_blast_run"
  )
})
