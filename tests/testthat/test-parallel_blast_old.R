testthat::test_that("new parallel_blast generates same output as old", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  withr::local_options(list(lifecycle_verbosity = "quiet"))

  ASVs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "asvs_test",
    ext = "txt"
  ))

  db_path <- tmp_blast_db_path

  blast_res_new <- parallel_blast(
    query_seqs = ASVs_test,
    db_path = db_path,
    total_cores = 1,
    perc_id = 80,
    num_threads = 1,
    perc_qcov_hsp = 80,
    num_alignments = 4,
    blast_type = "blastn",
    verbose = "silent"
  )

  # remove metadata
  attributes(blast_res_new)$BLASTr_metadata <- NULL

  blast_res_old <- parallel_blast_old(
    asvs = ASVs_test,
    db_path = db_path,
    out_file = NA,
    out_RDS = NA,
    total_cores = 1,
    perc_id = 80,
    num_threads = 1,
    perc_qcov_hsp = 80,
    num_alignments = 4,
    blast_type = "blastn",
    verbose = "silent"
  )

  testthat::expect_equal(
    dplyr::select(blast_res_new, -c("1_staxid")),
    dplyr::select(blast_res_old, -c("1_staxid"))
  )
  testthat::expect_s3_class(blast_res_new, "tbl_df")
})


testthat::test_that("parallel_blast_old - query with no hits", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  withr::local_options(list(lifecycle_verbosity = "quiet"))

  fail_query_seq <- "AGAGTACTACAAGTGCTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast_old(
    query_seqs = fail_query_seq,
    db_path = db_path,
    blast_type = "blastn"
  )

  testthat::expect_s3_class(blast_res, "tbl_df")
})

testthat::test_that("parallel_blast_old - query that errors", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  withr::local_options(list(lifecycle_verbosity = "quiet"))

  fail_query_seq <- ">XAVADADAA"

  db_path <- tmp_blast_db_path

  parallel_blast_old(
    query_seqs = fail_query_seq,
    db_path = db_path,
    blast_type = "blastn"
  ) |>
    testthat::expect_error()
})


testthat::test_that("OLD - parallel blast - single thread", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  withr::local_options(list(lifecycle_verbosity = "quiet"))

  asvs_test <- readLines(
    fs::path_package("BLASTr", "extdata", "asvs_test", ext = "txt")
  )

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast_old(
    asvs = asvs_test,
    db_path = db_path,
    num_threads = 1,
    total_cores = 1
  )

  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
  testthat::expect_equal(nrow(blast_res), length(asvs_test))
  testthat::expect_equal(nrow(blast_res), 6L)
  testthat::expect_equal(ncol(blast_res), 57L)
})

testthat::test_that("OLD - parallel blast - empty args", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  withr::local_options(list(lifecycle_verbosity = "warning"))

  testthat::expect_warning(
    blast_res <- parallel_blast_old(asvs = NULL, db_path = NULL)
  )

  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
})
