testthat::test_that("background job gives the same results as the foreground", {
  skip_if_no_blast()
  seqs_test <- readLines(fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt"))
  fg <- parallel_blast(seqs_test, tmp_blast_db_path, total_cores = 2, num_alignments = 4, verbose = "silent")

  job_dir <- withr::local_tempdir()
  fs::dir_delete(job_dir)
  job <- parallel_blast(
    seqs_test, tmp_blast_db_path,
    total_cores = 2, num_alignments = 4, verbose = "silent",
    background = TRUE, job_dir = job_dir
  )
  testthat::expect_s3_class(job, "blastr_job")
  bg <- blast_job_result(blast_job(job_dir), progress = FALSE)
  attr(bg, "BLASTr_metadata") <- NULL
  attr(fg, "BLASTr_metadata") <- NULL
  testthat::expect_equal(bg, fg)
  st <- blast_job_status(job)
  testthat::expect_true(st$finished)
  testthat::expect_equal(st$done, length(unique(seqs_test)))
})

testthat::test_that("job directory must be new", {
  skip_if_no_blast()
  d <- withr::local_tempdir()
  writeLines("x", fs::path(d, "other.txt"))
  parallel_blast("ACGT", tmp_blast_db_path, job_dir = d, background = TRUE, verbose = "silent") |>
    testthat::expect_error(class = "blastr_job_dir_exists_error")
  blast_job(d) |> testthat::expect_error(class = "blastr_not_a_job_error")
})

testthat::test_that("cancel and resume a background job", {
  skip_if_no_blast()
  seqs_test <- readLines(fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt"))
  job_dir <- fs::path(withr::local_tempdir(), "job")
  job <- parallel_blast(
    seqs_test, tmp_blast_db_path,
    batch_size = 1, total_cores = 1, num_alignments = 4, verbose = "silent",
    background = TRUE, job_dir = job_dir
  )
  blast_job_cancel(job)
  for (i in 1:100) {
    if (!BLASTr:::scheduler_alive(job) && blast_job_status(job)$finished) break
    Sys.sleep(0.2)
  }
  st <- blast_job_status(job)
  testthat::expect_true(st$finished)
  if (st$cancelled) {
    testthat::expect_lt(st$done - sum(BLASTr:::read_batches(job)$status == "cancelled"), st$total)
  }
  blast_job_resume(job)
  res <- blast_job_result(job, progress = FALSE)
  testthat::expect_equal(length(unique(res$Sequence)), length(unique(seqs_test)))
  testthat::expect_true(all(exit_codes(res)$exit_code == 0L))
})

testthat::test_that("searching several databases together", {
  skip_if_no_blast()
  seqs_test <- readLines(fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt"))
  res <- parallel_blast(
    seqs_test, c(tmp_blast_db_path, tmp_blast_db_path),
    num_alignments = 4, verbose = "silent"
  )
  testthat::expect_equal(length(unique(res$Sequence)), length(unique(seqs_test)))
  testthat::expect_equal(attr(res, "BLASTr_metadata")$db_path, c(tmp_blast_db_path, tmp_blast_db_path))
})
