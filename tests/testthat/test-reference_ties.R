# The tie rule by references (at_best_refs) and the species flag of
# blast_lca(ties = "at_best_refs"), on hits made by hand for queries of 200 bp.
ref_hits <- function(...) {
  r <- dplyr::bind_rows(lapply(list(...), tibble::as_tibble))
  r$identity <- 100 * r$nident / r$length
  r$qcovhsp <- round(100 * (r$`query end` - r$`query start` + 1) / 200)
  r$bitscore <- r$nident
  r$`subject header` <- r$subject
  r$blast_rank <- as.integer(stats::ave(seq_len(nrow(r)), r$query, FUN = seq_along))
  r
}
h <- function(query, subject, staxid, btop, length = 200, nident = 200, qstart = 1, qend = 200) {
  list(
    query = query, subject = subject, staxid = staxid, btop = btop, length = length, nident = nident,
    `query start` = qstart, `query end` = qend
  )
}
q200 <- strrep("A", 200)

testthat::test_that("btop_states: substitutions, gaps and insertions by query position", {
  st <- btop_states(c("10AG5", "5-A5", "5A-5", "16", NA), c(1, 1, 1, 1, 1))
  testthat::expect_equal(as.vector(st[[1]]), "G")
  testthat::expect_equal(names(st[[1]]), "p11")
  testthat::expect_equal(names(st[[2]]), "i6")
  testthat::expect_equal(as.vector(st[[2]]), "A")
  testthat::expect_equal(names(st[[3]]), "p6")
  testthat::expect_equal(as.vector(st[[3]]), "-")
  testthat::expect_length(st[[4]], 0L)
  testthat::expect_null(st[[5]])
  # differences between two hits on a region
  testthat::expect_equal(state_diff(st[[1]], st[[4]], 1, 16), 1L)
  testthat::expect_equal(state_diff(st[[1]], st[[1]], 1, 16), 0L)
  testthat::expect_equal(state_diff(st[[1]], st[[4]], 12, 16), 0L)
})

testthat::test_that("at_best_refs: references as close to the best one as the query, at least 0.5%", {
  hits <- ref_hits(
    # identical to the query
    h("q1", "a", "9606", "200"),
    # one difference: within the 0.5% floor of the identical one
    h("q1", "b", "9598", "100AG99", nident = 199),
    # two differences: beyond it
    h("q1", "c", "9580", "50AC49AT99", nident = 198),
    # the query differs from the best reference by one base (0.5%)
    h("q2", "d", "9606", "100AG99", nident = 199),
    # two bases from the query, one from the best reference: tied
    h("q2", "e", "9598", "100AG50AC48", nident = 198),
    # two bases from the query, three from the best reference: not tied
    h("q2", "f", "9580", "60AC40AG98", nident = 198),
    # covers less than half of the query: cannot be told apart, tied
    h("q2", "g", "9580", "90", length = 90, nident = 90, qend = 90)
  )
  long <- blast_long(hits, c("q1", "q2"), c(q200, q200), num_alignments = 10L)
  testthat::expect_equal(names(long), long_cols)
  r <- stats::setNames(long$at_best_refs, long$subject)
  testthat::expect_equal(r[c("a", "b", "c")], c(a = TRUE, b = TRUE, c = FALSE))
  testthat::expect_equal(r[c("d", "e", "f", "g")], c(d = TRUE, e = TRUE, f = FALSE, g = TRUE))
  # the default rule ties f (0.45 point of Similarity score from the best)
  testthat::expect_true(long$at_best[long$subject == "f"])
  testthat::expect_equal(long$btop[long$subject == "d"], "100AG99")
})

testthat::test_that("at_best_refs keeps the tied hits beyond num_alignments", {
  hits <- ref_hits(
    h("q1", "a", "9606", "200"),
    h("q1", "b", "9598", "30AC169", nident = 199),
    h("q1", "c", "9580", "100AG99", nident = 199)
  )
  # rank by bitscore puts c after b; keep one hit only
  long <- blast_long(hits, "q1", q200, num_alignments = 1L, rank_by = "bitscore", tolerance = 0)
  testthat::expect_setequal(long$subject, c("a", "b", "c"))
  testthat::expect_true(all(long$at_best_refs))
})

testthat::test_that("without btop, at_best_refs is NA and blast_lca() says why", {
  hits <- ref_hits(h("q1", "a", "9606", NA), h("q1", "b", "9598", NA, nident = 199))
  long <- blast_long(hits, "q1", q200)
  testthat::expect_true(all(is.na(long$at_best_refs)))
  testthat::expect_error(blast_lca(long, ties = "at_best_refs", taxonomy_db = lca_db()), class = "blastr_no_btop_error")
})

testthat::test_that("blast_lca(ties = 'at_best_refs'): ancestor of the references and species flag", {
  db <- lca_db()
  hits <- ref_hits(
    # A: the query 1% from its only Homo sapiens reference, Pan far away
    h("q1", "a1", "9606", "50AC49AT99", nident = 198),
    h("q1", "a2", "9598", "10AC10AC10AC10AC10AC10AC10AC10AC10AC10AC90", nident = 190),
    # B: the query 0.5% from Homo sapiens, whose references differ by 2%
    h("q2", "b1", "9606", "100AG99", nident = 199),
    h("q2", "b2", "9606", "100AG10AC10AC10AC10AC55", nident = 195),
    # C: the query 3% from Homo sapiens, whose references differ by 1%
    h("q3", "c1", "9606", "20AC20AC20AC20AC20AC20AC74", nident = 194),
    h("q3", "c2", "9606", "20AC20AC20AC20AC20AC20AC10AG10AG52", nident = 192)
  )
  # three different sequences of 200 bp
  long <- blast_long(hits, c("q1", "q2", "q3"), c(q200, sub("A$", "C", q200), sub("^A", "C", q200)), num_alignments = 10L)
  res <- blast_lca(long, ties = "at_best_refs", taxonomy_db = db)
  testthat::expect_equal(names(res)[11:12], c("Species flag", "Flagged species"))
  testthat::expect_equal(res$`LCA name`, rep("Homo sapiens", 3))
  testthat::expect_equal(res$`Species flag`, c("single reference", NA, "beyond spread"))
  testthat::expect_equal(res$`Flagged species`, c("Homo sapiens", NA, "Homo sapiens"))
  # b2 is not tied (2% from b1) but is kept: same taxid as the best reference
  testthat::expect_false(long$at_best_refs[long$subject == "b2"])
  # the default rule gives the same columns as before
  testthat::expect_equal(ncol(blast_lca(long, taxonomy_db = db)), 10L)
})
