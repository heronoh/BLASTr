hits_df <- function() {
  data.frame(
    query = c("s1", "s1", "s1", "s2"),
    blast_rank = c(1L, 2L, 3L, 1L),
    subject = c("A", "B", "C", "D"),
    `subject header` = c("a", "b", "c", "d"),
    staxid = c("1", "2", "1", "3"),
    identity = c(99, 100, 100, 97),
    qcovhsp = c(100, 100, 100, 90),
    bitscore = c(500, 480, 480, 300),
    check.names = FALSE
  )
}

testthat::test_that("blast_long ranks hits, marks ties and keeps sequences without hits", {
  h <- hits_df()
  h$`Similarity score` <- similarity_score(h$identity, h$qcovhsp)
  long <- blast_long(h, c("s1", "s2", "s3"), c("AAA", "CCC", "GGG"), num_alignments = 10L)
  testthat::expect_equal(names(long), long_cols)
  testthat::expect_equal(long$Sequence, c("AAA", "AAA", "AAA", "CCC", "GGG"))
  testthat::expect_equal(long$subject, c("B", "C", "A", "D", NA))
  testthat::expect_equal(long$rank, c(1L, 1L, 3L, 1L, NA))
  testthat::expect_equal(long$tie, c(TRUE, TRUE, FALSE, FALSE, NA))
  testthat::expect_equal(long$at_best, c(TRUE, TRUE, FALSE, TRUE, NA))
  testthat::expect_equal(long$blast_rank, c(2L, 3L, 1L, 1L, NA))
  testthat::expect_true(all(is.na(long$low_qcov)))
  long <- blast_long(h, c("s1", "s2", "s3"), c("AAA", "CCC", "GGG"), perc_qcov_hsp = 95)
  testthat::expect_equal(long$low_qcov, c(FALSE, FALSE, FALSE, TRUE, NA))
  testthat::expect_equal(blast_summary(long)$`Best low qcov`, c(FALSE, TRUE, NA))
})

testthat::test_that("blast_long keeps all hits at best beyond num_alignments", {
  h <- hits_df()
  long <- blast_long(h, c("s1", "s2"), c("AAA", "CCC"), num_alignments = 1L)
  testthat::expect_equal(long$subject, c("B", "C", "D"))
})

testthat::test_that("blast_summary, wide and as_blast_long agree", {
  h <- hits_df()
  long <- blast_long(h, c("s1", "s2", "s3"), c("AAA", "CCC", "GGG"))
  s <- blast_summary(long, search_alignments = 50L)
  testthat::expect_equal(s$Sequence, c("AAA", "CCC", "GGG"))
  testthat::expect_equal(s$`N hits`, c(3L, 1L, 0L))
  testthat::expect_equal(s$`Check ID`, c("check ID", "best ID", NA))
  testthat::expect_equal(s$`Taxids at best`, c("2;1", "3", NA))
  wide <- widen_blast_hits(long, s, 2L)
  testthat::expect_equal(wide$`1_subject`, c("B", "D", NA))
  testthat::expect_equal(wide$`1_blast_rank`, c(2L, 1L, NA))
  back <- as_blast_long(wide)
  testthat::expect_equal(back$subject, c("B", "C", "D", NA))
  testthat::expect_equal(back$blast_rank, c(2L, 3L, 1L, NA))
})

testthat::test_that("as_blast_long reads old wide results (k_indentity, BLAST order)", {
  wide <- data.frame(
    Sequence = c("AAA", "CCC"),
    `1_subject` = c("A", NA), `1_indentity` = c(99, NA), `1_qcovhsp` = c(100, NA), `1_staxid` = c("1", NA),
    `2_subject` = c("B", NA), `2_indentity` = c(100, NA), `2_qcovhsp` = c(100, NA), `2_staxid` = c("2", NA),
    check.names = FALSE
  )
  long <- as_blast_long(wide)
  testthat::expect_equal(long$subject, c("B", "A", NA))
  testthat::expect_equal(long$blast_rank, c(2L, 1L, NA))
  testthat::expect_equal(long$rank, c(1L, 2L, NA))
  testthat::expect_equal(blast_summary(long)$`N hits`, c(2L, 0L))
})

testthat::test_that("blast_cmd_line hides the sequences and quotes arguments", {
  cmd <- blast_cmd_line("blastn", c("-db", "/x/core_nt", "-outfmt", "6 qaccver saccver"), 32L)
  testthat::expect_equal(
    cmd,
    "blastn -db /x/core_nt -outfmt '6 qaccver saccver' -query <batch.fasta> -out <batch.tsv> -num_threads 32"
  )
})
