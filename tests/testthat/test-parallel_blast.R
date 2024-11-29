testthat::test_that("parallel blast single thread", {
  asvs_test <- readLines(fs::path_package("BLASTr", "extdata", "asvs_test.txt"))
  db_test <- fs::path_package("BLASTr", "extdata", "minimal_db_blast.fasta")

  blast_res <- parallel_blast(
    asvs = asvs_test,
    db_path = db_test
  )

  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
  testthat::expect_equal(nrow(blast_res), length(asvs_test))
  testthat::expect_equal(nrow(blast_res), 6L)
})
