test_that("Get single taxid", {
  testthat::skip_if_offline()
  res <- get_tax_by_taxID(organisms_taxIDs = "63221", verbose = FALSE)

  testthat::expect_s3_class(res, "tbl_df")

  testthat::expect_equal(res$query_taxID, "63221")

  testthat::expect_true(stringr::str_detect(res$Sci_name, "Homo sapiens"))

  testthat::expect_equal(ncol(res), 13L)
})

test_that("Get single taxid empty df", {
  testthat::skip_if_offline()
  res <- get_tax_by_taxID(organisms_taxIDs = "10000000", verbose = FALSE)

  testthat::expect_s3_class(res, "tbl_df")

  testthat::expect_equal(nrow(res), 0L)

  testthat::expect_equal(ncol(res), 13L)
})

test_that("Get parallel Single taxid", {
  testthat::skip_if_offline()
  res <- parallel_get_tax(organisms_taxIDs = "63221", verbose = FALSE)

  testthat::expect_equal(res$query_taxID, "63221")

  testthat::expect_true(stringr::str_detect(res$Sci_name, "Homo sapiens"))
})

test_that("Get parallel Multiple taxid", {
  testthat::skip_if_offline()
  res <- parallel_get_tax(
    organisms_taxIDs = c("63221", "10000000"),
    total_cores = 1,
    retry_times = 2,
    verbose = FALSE
  )

  testthat::expect_equal(res$query_taxID, "63221")

  testthat::expect_true(stringr::str_detect(res$Sci_name, "Homo sapiens"))
})


test_that("Get parallel single-thread Multiple taxid", {
  testthat::skip_if_offline()
  res <- parallel_get_tax(
    organisms_taxIDs = c("63221", "10000000"),
    total_cores = 2,
    retry_times = 2,
    verbose = FALSE
  )

  testthat::expect_equal(res$query_taxID, "63221")

  testthat::expect_true(stringr::str_detect(res$Sci_name, "Homo sapiens"))
})


test_that("`parallel_get_tax()` with `retry_times = 0` ", {
  testthat::skip_if_offline()
  res <- parallel_get_tax(
    organisms_taxIDs = c("63221", "10000000"),
    total_cores = 1,
    retry_times = 0,
    verbose = FALSE
  )

  testthat::expect_equal(res$query_taxID, "63221")

  testthat::expect_true(stringr::str_detect(res$Sci_name, "Homo sapiens"))
})
