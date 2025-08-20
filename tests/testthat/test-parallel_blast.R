testthat::test_that("parallel blast - single thread", {
  asvs_test <- readLines(
    fs::path_package("BLASTr", "extdata", "asvs_test", ext = "txt")
  )
  db_test <- fs::path_package(
    "BLASTr", "extdata", "minimal_db_blast",
    ext = "fasta"
  )

  blast_res <- parallel_blast(
    asvs = asvs_test,
    db_path = db_test,
    num_threads = 1,
    total_cores = 1
  )

  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
  testthat::expect_equal(nrow(blast_res), length(asvs_test))
  testthat::expect_equal(nrow(blast_res), 6L)
  testthat::expect_equal(ncol(blast_res), 57L)
})

testthat::test_that("parallel blast - empty args", {
  blast_res <- parallel_blast(asvs = NULL, db_path = NULL)
  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
})
