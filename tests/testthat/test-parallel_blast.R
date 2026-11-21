testthat::test_that("parallel_blast works", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  seqs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "seqs_test",
    ext = "txt"
  ))

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast(
    query_seqs = seqs_test,
    db_path = db_path,
    total_cores = 1,
    perc_id = 80,
    num_threads = 1,
    perc_qcov_hsp = 80,
    num_alignments = 4,
    retry_times = 0,
    blast_type = "blastn",
    output = "wide",
    verbose = "silent"
  )

  testthat::expect_s3_class(blast_res, "tbl_df")
  testthat::expect_equal(nrow(blast_res), length(seqs_test))
  testthat::expect_equal(nrow(blast_res), 6L)
  testthat::expect_equal(ncol(blast_res), 1L + 9L + 4L * 19L)

  exit_codes_df <- exit_codes(blast_res)
  testthat::expect_s3_class(exit_codes_df, "data.frame")
  testthat::expect_equal(dim(exit_codes_df), c(6L, 3L))
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
    retry_times = 0L,
    verbose = "silent"
  )
  testthat::expect_equal(nrow(blast_res), 1L)

  exit_codes_df <- exit_codes(blast_res)

  testthat::expect_s3_class(exit_codes_df, "data.frame")
  testthat::expect_equal(dim(exit_codes_df), c(1L, 3L))
  testthat::expect_equal(exit_codes_df$exit_code, 3L)
  testthat::expect_equal(exit_codes_df$query_seq, fail_query_seq)

  blast_res <- parallel_blast(
    query_seqs = fail_query_seq,
    db_path = db_path,
    blast_type = "blastn",
    retry_times = 1L,
    verbose = "silent"
  )
  testthat::expect_equal(nrow(blast_res), 1L)

  blast_res <- parallel_blast(
    query_seqs = fail_query_seq,
    db_path = db_path,
    blast_type = "blastn",
    retry_times = 10L,
    verbose = "silent"
  )
  testthat::expect_equal(nrow(blast_res), 1L)

  blast_cmd_res <- blast_cmd(
    query_str = fail_query_seq,
    db_path = db_path,
    verbose = "silent"
  )
  testthat::expect_s3_class(blast_res, "tbl_df")
})


testthat::test_that("parallel_blast - 2 threads works", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  seqs_test <- readLines(fs::path_package(
    "BLASTr",
    "extdata",
    "seqs_test",
    ext = "txt"
  ))

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast(
    query_seqs = seqs_test,
    db_path = db_path,
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


testthat::test_that("parallel_blast - single thread", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  seqs_test <- readLines(
    fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt")
  )

  db_path <- tmp_blast_db_path

  blast_res <- parallel_blast(
    query_seqs = seqs_test,
    db_path = db_path,
    num_threads = 1,
    total_cores = 1,
    num_alignments = 4,
    output = "wide",
    verbose = "silent"
  )

  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
  testthat::expect_equal(nrow(blast_res), length(seqs_test))
  testthat::expect_equal(nrow(blast_res), 6L)
  testthat::expect_equal(ncol(blast_res), 1L + 9L + 4L * 19L)
})

testthat::test_that("parallel blast - empty args", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()
  parallel_blast(query_seqs = NULL, db_path = NULL) |>
    expect_error(class = "blastr_no_query_seqs_error")

  parallel_blast(query_seqs = character(0), db_path = NULL) |>
    expect_error(class = "blastr_no_query_seqs_error")

  parallel_blast(query_seqs = "AAAA", db_path = NULL) |>
    expect_error(class = "blastr_no_db_path_error")
})

testthat::test_that("parallel_blast - depracated arguments", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  withr::local_options(list(lifecycle_verbosity = "warning"))

  seqs_test <- readLines(
    fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt")
  )

  db_path <- tmp_blast_db_path

  testthat::expect_warning(
    object = {
      blast_res <- parallel_blast(
        query_seqs = seqs_test,
        db_path = db_path,
        num_threads = 1,
        num_alignments = 4,
        output = "wide",
        verbose = "silent",
        out_file = withr::local_tempfile(fileext = ".csv")
      )
    }
  )

  testthat::expect_s3_class(blast_res, "data.frame")
  testthat::expect_s3_class(blast_res, "tbl_df")
  testthat::expect_equal(nrow(blast_res), length(seqs_test))
  testthat::expect_equal(nrow(blast_res), 6L)
  testthat::expect_equal(ncol(blast_res), 1L + 9L + 4L * 19L)
})

testthat::test_that("parallel_blast long output: one row per hit, best first", {
  testthat::skip_on_cran()
  seqs_test <- readLines(fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt"))
  res <- parallel_blast(seqs_test, tmp_blast_db_path, num_alignments = 4, verbose = "silent")
  testthat::expect_equal(names(res), long_cols)
  testthat::expect_setequal(unique(res$Sequence), unique(seqs_test))
  hit <- res[!is.na(res$rank), ]
  testthat::expect_true(all(hit$rank[!duplicated(hit$Sequence)] == 1L))
  testthat::expect_true(all(stats::ave(hit$rank, hit$Sequence, FUN = function(r) all(diff(r) >= 0)) == 1))
  s <- blast_summary(res)
  testthat::expect_equal(s$Sequence, unique(res$Sequence))
  testthat::expect_true("Taxids at best" %in% names(s))
  meta <- attr(res, "BLASTr_metadata")
  testthat::expect_match(meta$blast_cmd, "^blastn -db .* -query <batch.fasta>")
  testthat::expect_false(any(vapply(seqs_test, grepl, logical(1L), x = meta$blast_cmd, fixed = TRUE)))
  testthat::expect_type(meta$db_dates, "character")
  testthat::expect_false(is.na(meta$blast_version))
})
