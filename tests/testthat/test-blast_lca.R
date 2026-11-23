lca_db <- function() {
  dir <- tiny_taxdump(withr::local_tempdir(.local_envir = parent.frame()))
  make_taxonomy_db(file.path(withr::local_tempdir(.local_envir = parent.frame()), "tax.sqlite"), taxdump = dir) |>
    suppressMessages()
}
long_of <- function(seq, staxid, at_best = TRUE) {
  tibble::tibble(
    Sequence = seq, rank = 1L, tie = TRUE, at_best = at_best, low_qcov = FALSE,
    blast_rank = seq_along(seq), subject = paste0("s", seq_along(seq)), staxid = staxid
  )
}

testthat::test_that("blast_lca: strict and consensus ancestors", {
  db <- lca_db()
  # 9 Homo sapiens + 1 Pan troglodytes tied
  x <- long_of(rep("A", 10), c(rep("9606", 9), "9598"))
  strict <- blast_lca(x, consensus = 1, taxonomy_db = db)
  testthat::expect_equal(strict$`LCA name`, "Hominidae")
  testthat::expect_equal(strict$`LCA rank`, "family")
  testthat::expect_equal(strict$`LCA support`, 1)
  cons <- blast_lca(x, taxonomy_db = db)
  testthat::expect_equal(cons$`LCA taxid`, "9606")
  testthat::expect_equal(cons$`LCA name`, "Homo sapiens")
  testthat::expect_equal(cons$`LCA rank`, "species")
  testthat::expect_equal(cons$`LCA support`, 0.9)
  testthat::expect_equal(cons$`N tied`, 10L)
})

testthat::test_that("blast_lca: few tied hits need all of them", {
  db <- lca_db()
  x <- long_of(c("A", "A"), c("9606", "9598"))
  testthat::expect_equal(blast_lca(x, taxonomy_db = db)$`LCA name`, "Hominidae")
})

testthat::test_that("blast_lca: only hits at best vote; unknown taxids do not", {
  db <- lca_db()
  x <- long_of(c("A", "A", "A", "B", "C"), c("9606", "0", "9598", "0", NA), at_best = c(TRUE, TRUE, FALSE, TRUE, NA))
  res <- blast_lca(x, taxonomy_db = db)
  testthat::expect_equal(res$Sequence, c("A", "B", "C"))
  testthat::expect_equal(res$`LCA name`, c("Homo sapiens", NA, NA))
  testthat::expect_equal(res$`N tied`, c(2L, 1L, NA))
  testthat::expect_equal(res$`N tied no taxid`, c(1L, 1L, 0L))
})

testthat::test_that("blast_lca takes a ready taxonomy and checks consensus", {
  db <- lca_db()
  x <- long_of(rep("A", 3), c("9606", "9606", "9605"))
  tax <- get_tax(x, taxonomy_db = db)
  res <- blast_lca(x, named_first = FALSE, taxonomy = tax, taxonomy_db = "/does/not/exist")
  testthat::expect_equal(res$`LCA name`, "Homo")
  testthat::expect_equal(res$`LCA rank`, "genus")
  # a hit placed in the genus still votes
  testthat::expect_equal(blast_lca(x, taxonomy = tax)$`LCA name`, "Homo")
  # one named only above genus does not, when a species is tied with it
  y <- long_of(rep("A", 3), c("9606", "9606", "9604"))
  res <- blast_lca(y, taxonomy = get_tax(y, taxonomy_db = db))
  testthat::expect_equal(res$`LCA name`, "Homo sapiens")
  testthat::expect_equal(res$`N tied`, 3L)
  testthat::expect_equal(res$`N tied no taxid`, 0L)
  testthat::expect_equal(blast_lca(y, named_first = FALSE, taxonomy_db = db)$`LCA name`, "Hominidae")
  testthat::expect_error(blast_lca(x, consensus = 0.5, taxonomy = tax))
})

testthat::test_that("blast_lca: ancestor without rank and the ranked one above it", {
  db <- lca_db()
  x <- long_of(c("A", "A", "B", "C"), c("9606", "9580", "9606", "9598"))
  x$Sequence[4] <- "B"
  res <- blast_lca(x, consensus = 1, taxonomy_db = db)
  testthat::expect_equal(res$`LCA name`, c("Hominoidea", "Hominidae"))
  testthat::expect_equal(res$`LCA rank`, c(NA, "family"))
  testthat::expect_equal(res$`LCA ranked name`, c("Primates", "Hominidae"))
  testthat::expect_equal(res$`LCA ranked rank`, c("order", "family"))
  testthat::expect_equal(res$`LCA ranked taxid`, c("9443", "9604"))
  # a species is its own ranked ancestor
  sp <- blast_lca(long_of("A", "9606"), taxonomy_db = db)
  testthat::expect_equal(sp$`LCA ranked name`, "Homo sapiens")
  testthat::expect_equal(sp$`LCA ranked rank`, "species")
})

testthat::test_that("blast_lca: the rank columns of taxonomy set the ranked ancestor", {
  db <- lca_db()
  x <- long_of(c("A", "A"), c("9606", "9598"))
  tax <- get_tax(x, taxonomy_db = db)
  res <- blast_lca(x, consensus = 1, taxonomy = tax[, c("taxid", "rank", "lineage", "lineage_taxids", "genus", "order")])
  testthat::expect_equal(res$`LCA name`, "Hominidae")
  testthat::expect_equal(res$`LCA ranked name`, "Primates")
  testthat::expect_equal(res$`LCA ranked rank`, "order")
})

testthat::test_that("blast_lca: uninformative hits do not vote", {
  db <- lca_db()
  x <- long_of(rep("A", 3), c("9606", "9598", "9598"))
  x$informative <- c(TRUE, FALSE, FALSE)
  res <- blast_lca(x, taxonomy_db = db)
  testthat::expect_equal(res$`LCA name`, "Homo sapiens")
  testthat::expect_equal(res$`N tied`, 1L)
  # only open nomenclature tied: it votes
  y <- long_of(c("B", "B"), c("9605", "9596"))
  testthat::expect_equal(blast_lca(y, taxonomy_db = db)$`LCA name`, "Hominidae")
})

testthat::test_that("blast_lca: a species in open nomenclature is ranked at its genus", {
  dir <- tiny_taxdump(withr::local_tempdir())
  dmp <- function(...) paste0(paste(..., sep = "\t|\t"), "\t|")
  cat(dmp(1000001, 9605, "species"), "\n", sep = "", file = file.path(dir, "nodes.dmp"), append = TRUE)
  cat(dmp(1000001, "Homo sp. X1", "", "scientific name"), "\n", sep = "", file = file.path(dir, "names.dmp"), append = TRUE)
  db <- make_taxonomy_db(file.path(withr::local_tempdir(), "tax.sqlite"), taxdump = dir) |> suppressMessages()
  res <- blast_lca(long_of("A", "1000001"), taxonomy_db = db)
  testthat::expect_equal(res$`LCA name`, "Homo sp. X1")
  testthat::expect_equal(res$`LCA ranked name`, "Homo")
  testthat::expect_equal(res$`LCA ranked rank`, "genus")
  testthat::expect_true(grepl(open_nomenclature_pattern, "Chlamydiota bacterium hydrae JB", perl = TRUE))
  testthat::expect_false(grepl(open_nomenclature_pattern, "Mycobacterium tuberculosis", perl = TRUE))
  testthat::expect_false(grepl(open_nomenclature_pattern, "Picea abies", perl = TRUE))
})
