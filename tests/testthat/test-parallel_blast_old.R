testthat::test_that("new parallel_blast generates same output as old", {
  skip_if_no_blast()

  withr::local_options(list(lifecycle_verbosity = "quiet"))

  seqs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "seqs_test",
    ext = "txt"
  ))

  db_path <- tmp_blast_db_path

  blast_res_new <- parallel_blast(
    query_seqs = seqs_test,
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
    query_seqs = seqs_test,
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

  # same hits, whatever their order (the new results are ranked by BLASTr)
  key <- function(x) {
    x <- x[!is.na(x$subject), ]
    sort(paste(x$Sequence, x$subject, x$identity, x$qcovhsp))
  }
  testthat::expect_equal(key(blast_res_new), key(as_blast_long(blast_res_old)))
  testthat::expect_s3_class(blast_res_new, "tbl_df")
})


testthat::test_that("parallel_blast_old - query with no hits", {
  skip_if_no_blast()

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
  skip_if_no_blast()

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
  skip_if_no_blast()

  withr::local_options(list(lifecycle_verbosity = "quiet"))

  seqs_test <- readLines(
    fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt")
  )

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast_old(
    query_seqs = seqs_test,
    db_path = db_path,
    num_threads = 1,
    total_cores = 1
  )

  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
  testthat::expect_equal(nrow(blast_res), length(seqs_test))
  testthat::expect_equal(nrow(blast_res), 6L)
  testthat::expect_equal(ncol(blast_res), 57L)
})

testthat::test_that("OLD - parallel blast - empty args", {
  skip_if_no_blast()
  withr::local_options(list(lifecycle_verbosity = "warning"))

  testthat::expect_warning(
    blast_res <- parallel_blast_old(query_seqs = NULL, db_path = NULL)
  )

  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
})
