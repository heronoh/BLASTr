seqs_test <- function() {
  readLines(fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt"))
}
fasta_test <- function() {
  fs::path_package("BLASTr", "extdata", "minimal_db_blast", ext = "fasta")
}

testthat::test_that("search_alignments keeps the num_alignments best hits", {
  skip_if_no_blast()
  s <- seqs_test()
  r2 <- parallel_blast(s, tmp_blast_db_path, num_alignments = 2, verbose = "silent")
  r10 <- parallel_blast(s, tmp_blast_db_path, num_alignments = 10, verbose = "silent")
  meta <- attr(r2, "BLASTr_metadata")
  testthat::expect_equal(meta$search_alignments, 100L)
  testthat::expect_false(any(startsWith(names(r2), "3_")))
  for (k in 1:2) {
    testthat::expect_equal(r2[[paste0(k, "_subject")]], r10[[paste0(k, "_subject")]])
  }
})

testthat::test_that("task is passed to BLAST", {
  skip_if_no_blast()
  r <- parallel_blast(seqs_test(), tmp_blast_db_path, task = "dc-megablast", verbose = "silent")
  testthat::expect_equal(attr(r, "BLASTr_metadata")$task, "dc-megablast")
  testthat::expect_true(all(exit_codes(r)$exit_code == 0L))
})

testthat::test_that("parallel_blast_hsps returns every HSP with typed fields", {
  skip_if_no_blast()
  s <- seqs_test()
  h <- parallel_blast_hsps(s, tmp_blast_db_path, total_cores = 2, verbose = "silent")
  testthat::expect_s3_class(h, "tbl_df")
  testthat::expect_equal(names(h), c("Sequence", "hit", "hsp", BLASTr:::hsp_fields_default))
  testthat::expect_type(h$pident, "double")
  testthat::expect_type(h$nident, "double")
  testthat::expect_true(all(h$Sequence %in% s))
  testthat::expect_true(all(h$hit >= 1L & h$hsp >= 1L))
  # one first HSP per subject; hits numbered 1, 2, ... in each sequence
  testthat::expect_false(any(duplicated(h[c("Sequence", "hit", "hsp")])))
  ranks <- tapply(h$hit, h$Sequence, function(x) identical(sort(unique(x)), seq_len(max(x))))
  testthat::expect_true(all(ranks))

  # no filters by default: at least the hits parallel_blast() finds
  w <- parallel_blast(s, tmp_blast_db_path, verbose = "silent")
  best <- h[h$hit == 1L & h$hsp == 1L, ]
  k <- match(unique(w$Sequence[!is.na(w$subject)]), best$Sequence)
  testthat::expect_false(anyNA(k))

  # chosen fields and filters
  h2 <- parallel_blast_hsps(s, tmp_blast_db_path, fields = c("saccver", "pident", "sstrand"),
                            perc_id = 95, max_hsps = 1, verbose = "silent")
  testthat::expect_equal(names(h2), c("Sequence", "hit", "hsp", "saccver", "pident", "sstrand"))
  testthat::expect_true(all(h2$pident >= 95))
  testthat::expect_true(all(h2$hsp == 1L))
})

testthat::test_that("parallel_blast_hsps searches a FASTA with -subject", {
  skip_if_no_blast()
  withr::local_options(BLASTr.log_dir = NULL)
  withr::local_envvar(BLASTR_LOG_DIR = "")
  log_dir <- withr::local_tempdir()
  s <- seqs_test()
  h <- parallel_blast_hsps(s, subject = fasta_test(), total_cores = 2, task = "blastn",
                           extra_args = c("-dust", "no"), verbose = "silent", log_dir = log_dir)
  testthat::expect_gt(nrow(h), 0L)
  meta <- attr(h, "BLASTr_metadata")
  testthat::expect_true(is.na(meta$strategy$threads))
  # no "num_threads is ignored" warning: -num_threads is not passed
  testthat::expect_false(any(grepl("num_threads", exit_codes(h)$stderr)))
  lg <- read_blast_log(log_dir)
  testthat::expect_equal(lg$`function`, "parallel_blast_hsps")
  testthat::expect_equal(lg$task, "blastn")
  testthat::expect_true(is.na(lg$pct_best_id))
})

testthat::test_that("parallel_blast_hsps checks its target", {
  parallel_blast_hsps("ACGT") |> testthat::expect_error(class = "blastr_no_db_path_error")
  parallel_blast_hsps("ACGT", db_path = "x", subject = "y") |>
    testthat::expect_error(class = "blastr_no_db_path_error")
  parallel_blast_hsps("ACGT", subject = tempfile()) |>
    testthat::expect_error(class = "blastr_no_db_path_error")
})

testthat::test_that("log rewritten with the union of columns when they change", {
  skip_if_no_blast()
  withr::local_options(BLASTr.log_dir = NULL)
  withr::local_envvar(BLASTR_LOG_DIR = "")
  log_dir <- withr::local_tempdir()
  path <- blast_log_path(log_dir)
  fs::dir_create(fs::path_dir(path))
  writeLines(c("run_id\tn_seqs\told_col", "r0\t3\tx"), path)
  parallel_blast(seqs_test(), tmp_blast_db_path, verbose = "silent", log_dir = log_dir)
  lg <- read_blast_log(log_dir)
  testthat::expect_equal(nrow(lg), 2L)
  testthat::expect_true(all(c("old_col", "function", "search_alignments") %in% names(lg)))
  testthat::expect_equal(lg$run_id[1], "r0")
})
