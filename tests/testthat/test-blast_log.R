testthat::test_that("no log unless a log directory is set", {
  withr::local_options(BLASTr.log_dir = NULL)
  withr::local_envvar(BLASTR_LOG_DIR = "")
  testthat::expect_null(blast_log_path())
  withr::local_envvar(BLASTR_LOG_DIR = "/tmp/x")
  testthat::expect_equal(as.character(blast_log_path()), "/tmp/x/BLASTr_logs/blastr_runs.tsv")
  withr::local_options(BLASTr.log_dir = "/tmp/y")
  testthat::expect_equal(as.character(blast_log_path()), "/tmp/y/BLASTr_logs/blastr_runs.tsv")
  testthat::expect_equal(as.character(blast_log_path("/tmp/z")), "/tmp/z/BLASTr_logs/blastr_runs.tsv")
})

testthat::test_that("foreground and background searches add one line each", {
  testthat::skip_on_cran()
  withr::local_options(BLASTr.log_dir = NULL)
  withr::local_envvar(BLASTR_LOG_DIR = "")
  seqs_test <- readLines(fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt"))
  log_dir <- withr::local_tempdir()

  res <- parallel_blast(
    c(seqs_test, seqs_test[1]), tmp_blast_db_path,
    total_cores = 2, verbose = "silent", log_dir = log_dir
  )
  meta <- attr(res, "BLASTr_metadata")
  testthat::expect_true(fs::file_exists(meta$log_file))
  testthat::expect_gte(meta$cpu_s, 0)

  job_dir <- fs::path(withr::local_tempdir(), "job")
  job <- parallel_blast(
    seqs_test, tmp_blast_db_path,
    total_cores = 2, verbose = "silent", log_dir = log_dir,
    background = TRUE, job_dir = job_dir
  )
  blast_job_result(job, progress = FALSE)
  blast_job_result(job, progress = FALSE) # collected twice, logged once

  lg <- read_blast_log(log_dir)
  testthat::expect_equal(nrow(lg), 2L)
  testthat::expect_equal(lg$mode, c("foreground", "background"))
  testthat::expect_equal(lg$n_input, c(length(seqs_test) + 1L, length(seqs_test)))
  testthat::expect_equal(lg$n_seqs, rep(length(unique(seqs_test)), 2L))
  testthat::expect_equal(lg$len_max, rep(max(nchar(unique(seqs_test))), 2L))
  testthat::expect_true(all(lg$cpu_s >= 0))
  testthat::expect_true(all(lg$pct_with_hits >= 0 & lg$pct_with_hits <= 100))
  testthat::expect_equal(lg$pct_best_id[1], lg$pct_best_id[2])
  testthat::expect_false(anyNA(lg$blast_version))
})
