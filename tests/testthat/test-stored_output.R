# A real BLAST+ output (inst/extdata/blast_output_example.tsv, written by
# data-raw/blast_output_example.R) checks the steps after the search where
# BLAST+ is not installed.
stored_search <- function() {
  seqs <- unique(readLines(fs::path_package("BLASTr", "extdata", "seqs_test", ext = "txt")))
  list(
    seqs = seqs,
    ids = sprintf("s%06d", seq_along(seqs)),
    hits = read_blast_hits(fs::path_package("BLASTr", "extdata", "blast_output_example", ext = "tsv"))
  )
}

testthat::test_that("a stored BLAST+ output is read, ranked and summarised without BLAST+", {
  x <- stored_search()
  testthat::expect_equal(nrow(x$hits), 16L)
  testthat::expect_equal(as.vector(table(x$hits$query)), c(1L, 1L, 1L, 1L, 12L))
  long <- blast_long(x$hits, x$ids, x$seqs, num_alignments = 10L,
                     uninformative = uninformative_references(), perc_qcov_hsp = 80)
  testthat::expect_equal(names(long), long_cols)
  testthat::expect_equal(unique(long$Sequence), x$seqs)
  # the second sequence has no hits: one row with NA hit fields
  testthat::expect_true(is.na(long$subject[long$Sequence == x$seqs[2]]))
  # the sixth: 10 ranked hits kept, the environmental one after them, unranked
  s6 <- long[long$Sequence == x$seqs[6], ]
  testthat::expect_equal(sum(!is.na(s6$rank)), 10L)
  testthat::expect_equal(s6$subject[1], "KX381515.1")
  testthat::expect_equal(s6$uninformative_category[is.na(s6$rank)], "environmental")
  testthat::expect_equal(s6$rank[s6$subject %in% c("KX381509.1", "KX381634.1")], c(4L, 4L))
  # summary: without the counts of the search, from the rows; with all the hits, 12
  s <- blast_summary(long, search_alignments = 100L)
  testthat::expect_equal(s$`N hits`, c(1L, 0L, 1L, 1L, 1L, 11L))
  testthat::expect_equal(s$`N uninformative`, c(0L, 0L, 0L, 0L, 0L, 1L))
  testthat::expect_equal(s$`Check ID`, c("check ID", NA, "check ID", "check ID", "check ID", "check ID"))
  full <- blast_summary(blast_long(x$hits, x$ids, x$seqs, uninformative = uninformative_references()))
  testthat::expect_equal(full$`N hits`[6], 12L)
})

testthat::test_that("the stored output is what BLAST+ gives now", {
  skip_if_no_blast()
  x <- stored_search()
  res <- parallel_blast(x$seqs, tmp_blast_db_path, verbose = "silent")
  long <- blast_long(x$hits, x$ids, x$seqs, num_alignments = 10L,
                     uninformative = uninformative_references(), perc_qcov_hsp = 80)
  testthat::expect_equal(res$subject, long$subject)
  testthat::expect_equal(res$identity, long$identity)
  testthat::expect_equal(res$rank, long$rank)
})
