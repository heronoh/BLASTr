testthat::test_that("new parallel_blast generates same output as old", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
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
    dplyr::select(blast_res_new, -c("1_staxid", "exit_code", "stderr")),
    dplyr::select(blast_res_old, -c("1_staxid"))
  )
  testthat::expect_s3_class(blast_res_new, "tbl_df")
})


testthat::test_that("parallel_blast - query with no hits", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  fail_query_seq <- "AGAGTACTACAAGTGCTAGCCATAAACTTAAATGAAGCTATACTAAACTCGTTCGCCAAAGCTTAAAACTCATAGGACTTGGCGGTGTTTCAGACCCAC"

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast(
    query_seqs = fail_query_seq,
    db_path = db_path,
    blast_type = "blastn"
  )

  blast_cmd_res <- blast_cmd(
    query_str = fail_query_seq,
    db_path = db_path
  )
  testthat::expect_s3_class(blast_res, "tbl_df")
})

testthat::test_that("parallel_blast - query that errors", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  fail_query_seq <- ">XAVADADAA"

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast(
    query_seqs = fail_query_seq,
    db_path = db_path,
    blast_type = "blastn",
    retry_times = 0L
  )
  testthat::expect_equal(nrow(blast_res), 1L)

  blast_res <- parallel_blast(
    query_seqs = fail_query_seq,
    db_path = db_path,
    blast_type = "blastn",
    retry_times = 1L
  )
  testthat::expect_equal(nrow(blast_res), 2L)

  blast_res <- parallel_blast(
    query_seqs = fail_query_seq,
    db_path = db_path,
    blast_type = "blastn",
    retry_times = 10L
  )
  testthat::expect_equal(nrow(blast_res), 11L)

  blast_cmd_res <- blast_cmd(
    query_str = fail_query_seq,
    db_path = db_path
  )
  testthat::expect_s3_class(blast_res, "tbl_df")
})


testthat::test_that("OLD - parallel_blast works", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  ASVs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "asvs_test",
    ext = "txt"
  ))

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast_old(
    asvs = ASVs_test,
    db_path = db_path,
    out_file = NA,
    out_RDS = NA,
    total_cores = 2,
    perc_id = 80,
    num_threads = 1,
    perc_qcov_hsp = 80,
    num_alignments = 4,
    blast_type = "blastn",
    verbose = "silent"
  )

  testthat::expect_s3_class(blast_res, "tbl_df")
})


testthat::test_that("OLD - parallel blast - single thread", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
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
  testthat::expect_warning(
    blast_res <- parallel_blast_old(asvs = NULL, db_path = NULL)
  )
  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
})
