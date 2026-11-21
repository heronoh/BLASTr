testthat::test_that("similarity_score weights identity 10x", {
  testthat::expect_equal(similarity_score(100, 100), 100)
  testthat::expect_equal(similarity_score(98, 87), (980 + 87) / 11)
  testthat::expect_equal(similarity_score(c(90, NA), c(80, 80)), c((900 + 80) / 11, NA))
})

hits_of <- function(identity, staxid, query = "q1") {
  data.frame(query = query, rank = seq_along(identity), identity = identity, staxid = staxid)
}

testthat::test_that("Check ID: same taxid in the top hits is best ID", {
  res <- check_hits(hits_of(c(100, 100, 99.8), c("1", "1", "1")))
  testthat::expect_equal(res$`Check ID`, "best ID")
  testthat::expect_equal(res$`N hits`, 3L)
  testthat::expect_equal(res$`N taxids`, 1L)
})

testthat::test_that("Check ID: other taxid tied in identity is check ID", {
  res <- check_hits(hits_of(c(100, 99.7, 95), c("1", "2", "1")))
  testthat::expect_equal(res$`Check ID`, "check ID")
  testthat::expect_equal(res$`N taxids at best identity`, 2L)
})

testthat::test_that("Check ID: clearly less identical hits are another taxon", {
  # 100% vs 99%: the second hit is most likely another species
  res <- check_hits(hits_of(c(100, 99, 98), c("1", "2", "3")))
  testthat::expect_equal(res$`Check ID`, "best ID")
  testthat::expect_equal(res$`N taxids`, 3L)
  testthat::expect_equal(res$`N hits same taxid as 1st`, 1L)
  testthat::expect_equal(res$`N taxids at best identity`, 1L)
})

testthat::test_that("Check ID: first hit less identical than another is check ID", {
  # BLAST ranks by score; a shorter alignment can be ranked lower but be more identical
  res <- check_hits(hits_of(c(99, 100), c("1", "2")))
  testthat::expect_equal(res$`Check ID`, "check ID")
})

testthat::test_that("Check ID: only the first check_n hits count", {
  h <- hits_of(c(100, 100, 100, 100), c("1", "1", "1", "2"))
  testthat::expect_equal(check_hits(h)$`Check ID`, "best ID")
  testthat::expect_equal(check_hits(h, check_n = 4)$`Check ID`, "check ID")
  testthat::expect_equal(check_hits(h)$`N taxids at best identity`, 2L)
})

testthat::test_that("Check ID: tolerance is configurable", {
  h <- hits_of(c(100, 99.7), c("1", "2"))
  testthat::expect_equal(check_hits(h, tolerance = 0.5)$`Check ID`, "check ID")
  testthat::expect_equal(check_hits(h, tolerance = 0.2)$`Check ID`, "best ID")
})

testthat::test_that("check_hits handles several queries and no hits", {
  h <- rbind(hits_of(c(100, 100), c("1", "2"), "a"), hits_of(100, "3", "b"))
  res <- check_hits(h)
  testthat::expect_equal(res$query, c("a", "b"))
  testthat::expect_equal(res$`Check ID`, c("check ID", "best ID"))
  testthat::expect_equal(nrow(check_hits(h[0, ])), 0L)
})

testthat::test_that("Check ID: unknown taxids (0) do not confirm agreement", {
  # databases built without a taxid map report staxid 0
  res <- check_hits(hits_of(c(100, 100, 99), c("0", "0", "0")))
  testthat::expect_equal(res$`Check ID`, "check ID")
  testthat::expect_equal(res$`N taxids`, 0L)
  testthat::expect_true(is.na(res$`N hits same taxid as 1st`))
  res <- check_hits(hits_of(c(100, 95), c("1", "0")))
  testthat::expect_equal(res$`Check ID`, "best ID")
  testthat::expect_equal(res$`N taxids`, 1L)
})
