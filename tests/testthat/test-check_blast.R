test_that("Check if BLAST CLI is available", {
  expect_equal(class(check_bin("blastn")), "character")
})
